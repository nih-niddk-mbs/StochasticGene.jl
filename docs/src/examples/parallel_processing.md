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
