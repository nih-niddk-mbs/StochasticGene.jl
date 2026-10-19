# Parallel And Batch Processing

This page summarizes current options for running many fits. For complete
cluster examples, see [Cluster and batch workflows](../cluster_batch_workflows.md).

## Parallel Chains In One Fit

For Metropolis-Hastings, `nchains` controls how many chains are run and merged.
When using distributed workers, start Julia with `-p` to match the number of
parallel chains:

```bash
julia -t 1 -p 4 fitscript_rna-HCT116_MOCK_2001.jl CENPL
```

The corresponding fit call should use:

```julia
fit(;
    nchains = 4,
    parallel = :distributed,
    datatype = "rna",
    datapath = "HCT116_testdata",
    gene = "CENPL",
    datacond = "MOCK",
    resultfolder = "HCT116_test",
)
```

## Batched Coupled Trace Likelihood

For supported hierarchical coupled trace models sharing a transition matrix,
`fit(...; likelihood_executor=:batched)` batches the forward recursion over
traces. This is separate from distributed chains and from Swarm job bundling.
The default executor remains `:fast`.

With several chains on separate workers, use one BLAS thread per worker.
For example, when continuing a keyed fit with its saved info and rate files:

```julia
@everywhere using StochasticGene, LinearAlgebra
@everywhere LinearAlgebra.BLAS.set_num_threads(1)

fit(;
    key = "H2-targets135-free-h",
    root = "/home/carsonc/nadya",
    resultfolder = "3Prime-hidden-multitarget-2026-10-16-restart",
    nchains = 16,
    parallel = :distributed,
    inference_method = :mh,
    likelihood_executor = :batched,
    maxtime = 1800.0,
)
```

Start the script with `julia -t 1 -p 16 fitscript.jl` (and `--project=...` if
using a checkout). Each worker runs one chain with one BLAS thread. On Biowulf,
request CPUs for all chains, as for your existing distributed fits.

This executor batches 16 traces within each chain; it does not require 16 Julia
threads and is unrelated to Swarm's `-b` job bundling. It reuses the shared
transition matrix in a batched forward recursion and stores only two forward
frames at a time. The cache described below is available in the current
checkout; released 2.0.7 rebuilds emissions at each likelihood call.

The default executor remains `:fast`, which can also be selected explicitly to
return to the reference implementation. Scalar/batched likelihood agreement
and short distributed fits have been checked on H2/H3-related test cases;
speed depends on model size and hardware, and this is not a convergence
guarantee. The regression tests live in `test/batched_forward_tests.jl`; the
read-only real-data benchmark and temporary-output fit check are
`test/benchmark_hidden_likelihood.jl` and `test/batched_fit_smoke.jl`.
The batched MH path reuses forward buffers within each chain. For coupled
reporters using `prob_Gaussian`, it also caches emission probabilities for
distinct reporter-count combinations rather than every joint state. Changes
to noise parameters, observations, or reporter definitions invalidate the
affected entries. Custom emission functions are evaluated without this
emission cache. Each chain owns its cache, released at the end of the chain;
the transition solver and statistical model are unchanged.

For repeated direct likelihood calls, the same cache can be scoped explicitly:

```julia
StochasticGene.with_batched_hmm_cache() do
    loglikelihood(param, data, model; hmm_stack=:batched)
    loglikelihood(next_param, data, model; hmm_stack=:batched)
end
```

Without this scope, direct calls still work but do not retain a cache across
evaluations. Performance depends on the model and hardware; see
`test/benchmark_batched_cache.jl` for a comparison against uncached batching.

### Optional Startup Tuning

Add these keywords to an existing MH `fit` call:

```julia
likelihood_executor = :batched,
hmm_batchsize = 16,
tune_hmm_batchsize = true,
hmm_batch_tuning_seconds = 120.0,
```

Each chain tunes independently at its initialized parameter vector, even with
`warmupsteps=0`. It warms candidate sizes 8, 16, 32, 64 and the requested initial
size, then collects up to three balanced timing rounds in alternating order.
At least two complete rounds are required. The smallest size within 3% of the
fastest median is selected and frozen before sampler warmup/sampling.
The chosen size, elapsed time, completed rounds, and status are printed in the
worker log. Requested settings are retained in the run spec and generated
fitscripts; the per-chain selected sizes are reported in logs, not the run spec.

Tuning is capped at the smaller of `hmm_batch_tuning_seconds` and 10% of
`maxtime`; its elapsed time is deducted from the remaining MH budget. This is a
soft deadline because an in-flight likelihood evaluation cannot be interrupted.
With too little time, the original size is kept. Unsupported paths skip tuning.
Likelihood disagreement aborts the fit rather than silently choosing a size.
No sampler random numbers are used, and no proposal adaptation is introduced.

The default is `tune_hmm_batchsize=false`. Set `hmm_batchsize` explicitly to
reuse a known choice without tuning. Tuning requires `likelihood_executor=:batched`.
Workers tune concurrently without a synchronization barrier, so noisy timing
under contention can select different sizes across chains. Benchmark on the
intended hardware; a Mac result is not a Biowulf performance guarantee.

### Standalone Experiment

The repository also provides a read-only batch-size experiment:

```bash
julia --project=. test/tune_hidden_batchsize.jl /path/to/nadya H2-targets135-free-h H3-targets135-free-h
```

This script reads the hidden-model fits from
`results/3Prime-hidden-multitarget-2026-10-16-restart` under the supplied root.
It warms sizes 8, 16, 32, and 64, compares three timings per size in randomized
order, and selects the smallest size within 3% of the fastest median. It then
compares the selection with size 16 on three fresh parameter proposals,
checking both total and per-trace likelihood agreement. It runs one process
with one BLAS thread and does not sample or modify source results.

This is computational tuning, not MCMC warmup or proposal adaptation. Default
fits still use size 16 and do not tune unless explicitly enabled. A selection
on one machine is not a recommendation for a different CPU or for concurrent
chains; benchmark under the intended workload before changing production runs.

## Optional Joint Residence Regularization

On the experimental `feature/joint-residence-regularizer` branch, coupled MH
trace fits can use an additional stationary ON/OFF occupancy penalty. The
default `regularizers=()` leaves the objective unchanged. This currently
requires one shared stationary distribution (including hierarchical models
with shared kinetics); grid, shared-dataset, and varying-kinetics paths are
not supported. NUTS, ADVI and `optimize=true` reject active regularizers.

```julia
# Rows/columns: ON_ON, ON_OFF, OFF_ON, OFF_OFF; enhancer is unit 2 here.
target = [0.041509, 0.210105, 0.108163, 0.640223]
# Example scale only: replace with the paired-trace bootstrap covariance.
covariance = 0.001 .* (Matrix{Float64}(I, 4, 4) .- 0.25)
occupancy = JointResidenceRegularizer(target;
    covariance=covariance, strength=1.0, units=(2, 1))

# Add to the existing fully specified fit call:
# regularizers=(occupancy,)
```

The four fractions sum to one, so they are not independent. The penalty uses
a covariance pseudoinverse, computed once per chain. ON means reporter count
greater than zero. State masks are also cached, and the HMM supplies its
already-computed stationary vector; no second stationary solve is performed
for each proposal. Regularizers and caches are scoped separately per chain.

MH samples `loglikelihood + logprior - penalty` (subject to the existing
temperature). The sampled `ll`, ML selection, AIC inputs and pointwise WAIC
inputs remain unpenalized. WAIC is nevertheless computed under a regularized
sampling distribution, not the original posterior. `regularization_<key>.csv`
reports fractions, individual/total penalties and the composite objective at
the raw-ML and coordinate-median points. Settings are retained in the run info
and fully generated fitscripts.

Use `regularization_report(param, data, model; regularizers=(occupancy,))`
for a read-only evaluation before fitting. Zero strength disables a term.
Targets derived from the same traces deliberately reuse information; they
are not independent FISH data. Stationary occupancy is not trace-centered
covariance. Hard-Viterbi target fractions may be biased relative to latent
occupancy; validate strengths using simulated decoding and held-out traces.

## Many Genes

For scRNA/RNA histogram sweeps, generate one shared fitscript and one command
file with one command per gene:

```julia
out = makeswarm_genes(
    root = "/data/carsonc/scrna",
    datapath = "RamosNELFA_NEG_IFNa_rep1_Sdata",
    datacond = "NIr1",
    cell = "U3A",
    resultfolder = "RamosNELFA_NEG_IFNa_rep1",
    filedir = "run-RamosNELFA_NEG_IFNa_rep1",
    nchains = 4,
    nthreads = 1,
)
```

By default, folder scanning uses `checkgenes` and keeps genes with available
halflife/allele metadata. Pass `filter_metadata=false` only when you want to
scan filenames directly.

## Biowulf Swarm

On Biowulf, submit the single command file with `swarm`. For large gene panels,
use `-b` to bundle several gene commands into each Slurm array task:

```bash
cd run-RamosNELFA_NEG_IFNa_rep1
swarm -f fit.swarm -b 20 -g 24 -t 4 --time 25:00 --merge-output --module julia
```

Here each Julia command gets `-t 4` CPUs and 24 GB. The `-b 20` flag means each
array task runs 20 gene commands serially; it reduces scheduler overhead without
changing per-command resources.

## Slurm Or GNU Parallel

On systems without Biowulf `swarm`, ask `makeswarm_genes` to write a launcher:

```julia
makeswarm_genes(
    datapath = "HCT116_testdata",
    datacond = "MOCK",
    resultfolder = "HCT116_test",
    filedir = "run-HCT116-testdata-rna",
    scheduler = :slurm,
    scheduler_jobs = 100,
    slurm_mem = "8G",
    slurm_time = "02:00:00",
)
```

Submit with:

```bash
sbatch run-HCT116-testdata-rna/fit_slurm.sh
```

For GNU Parallel:

```julia
makeswarm_genes(
    datapath = "HCT116_testdata",
    datacond = "MOCK",
    resultfolder = "HCT116_test",
    filedir = "run-HCT116-testdata-rna",
    scheduler = :parallel,
    scheduler_jobs = 16,
)
```

Run with:

```bash
./run-HCT116-testdata-rna/fit_parallel.sh
```

## Sequential Local Fallback

The generated `fit.swarm` file is a plain shell command list, so a regular
terminal can run it sequentially:

```bash
bash fit.swarm
```

This is useful for smoke tests or small panels, but it will not parallelize
across genes.
