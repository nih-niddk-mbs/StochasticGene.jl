# RNA Histogram Analysis

This example shows the current `datatype="rna"` workflow for steady-state RNA
count histograms from smFISH or scRNA-seq.

## Data Layout

For `datatype="rna"`, StochasticGene expects one text file per gene/condition:

```text
data/HCT116_testdata/
├── CENPL_MOCK.txt
└── MYC_MOCK.txt
```

The filename pattern is:

```text
GENE_COND.txt
```

The first column is the histogram over RNA copy number bins. For example,
`MYC_MOCK.txt` is loaded with:

```julia
gene = "MYC"
datacond = "MOCK"
datapath = "data/HCT116_testdata"
```

Paths are resolved the same way as `fit`: if needed, a relative `datapath` is
looked up under `root/data`.

Choose the histogram support policy explicitly when comparing analyses:

- `rna_truncation=:legacy` (default) reproduces v0.7.8 by retaining 99% of
  counts, capped at 1000 bins.
- `rna_truncation=:none` retains the complete histogram.

The setting applies to direct RNA fits, RNA legs of combined data, and
`makeswarm_genes` scripts.

## Single-Gene Fit

```julia
using StochasticGene

fits, stats, measures, data, model, options = fit(
    G = 2,
    R = 0,
    S = 0,
    insertstep = 1,
    transitions = ([1, 2], [2, 1]),
    datatype = "rna",
    datapath = "data/HCT116_testdata",
    gene = "MYC",
    cell = "HCT116",
    datacond = "MOCK",
    resultfolder = "HCT116_test",
    nchains = 2,
    maxtime = 60.0,
    samplesteps = 1_000_000,
    propcv = 0.01,
    rna_truncation = :legacy,
    burst = true,
    optimize = true,
)
```

If `decayrate=-1.0` (the default), `fit` tries to look up gene-specific decay
information from `data/halflives` under `root`, converts the halflife to a
decay rate in inverse minutes, and uses `1.0` if no value can be resolved.
Negative decay placeholders in an explicit `priormean` are replaced before
the starting vector is constructed, even when another run supplies starting
rates. If allele metadata exist under `data/alleles`, `fit` can also use those
to set `nalleles`.

`burst=true` and `optimize=true` request additional fit-time analyses:

- `burst=true` computes burst size over posterior samples and writes
  `burst_*.txt`.
- `optimize=true` starts LBFGS from the best sampled point and writes
  `optimized_*.txt`.

These files cannot be reconstructed by `make_dataframes` if they were not
created during fitting.

## Interpreting Rate, Burst, And Optimizer Outputs

Each `rates_*.txt` contains a header followed by four complete model-rate rows:

1. sampled ML: the largest likelihood encountered during inference,
2. posterior mean,
3. posterior median,
4. last stored posterior sample.

The sampled-ML row is not the output of `optimize=true`. Numerical optimization
has its own file and assembled columns. `OptimizerConverged=true` means the
optimizer's stopping test passed; it does not prove a unique solution. Models
with flat likelihood ridges can have nearly identical optimized objectives but
very different individual optimized rates.

For a G=2 RNA model, posterior burst size is evaluated sample by sample as:

```math
\mathrm{burst\ size} = \frac{\mathrm{Eject}}{\mathrm{Rate21}}.
```

The output reports `BurstMean`, `BurstSD`, `BurstMedian`, and `BurstMAD`
(the raw file also stores quantiles). Because the ratio is evaluated within
each posterior sample, `BurstMedian` is not generally
`median(Eject) / median(Rate21)`. Derived burst size can also be much more
stable than either rate separately when the posterior follows a correlated
likelihood ridge.

## Genome-Scale scRNA Sweeps

For many genes, use `makeswarm_genes`. The folder-scanning form infers genes
from filenames matching `GENE_COND.txt`, writes one shared fitscript using
`ARGS[1]` as the gene name, and writes one command per gene.

```julia
out = makeswarm_genes(
    root = "/data/carsonc/scrna",
    datapath = "RamosNELFA_NEG_IFNa_rep1_Sdata",
    datacond = "NIr1",
    cell = "U3A",
    resultfolder = "RamosNELFA_NEG_IFNa_rep1",
    filedir = "run-RamosNELFA_NEG_IFNa_rep1",
    G = 2,
    R = 0,
    S = 0,
    insertstep = 1,
    transitions = ([1, 2], [2, 1]),
    nchains = 2,
    nthreads = 1,
    maxtime = 3600.0,
    samplesteps = 1_000_000,
    propcv = 0.01,
    rna_truncation = :legacy,
    burst = true,
    optimize = true,
    project = "/home/carsonc/github/StochasticGene.jl/",
)
```

By default, folder scanning uses `checkgenes`, so the gene list is restricted to
genes with the halflife/allele metadata needed by the RNA fit. Pass
`filter_metadata=false` only when you deliberately want to scan filenames
without that metadata gate.

`makeswarm_genes` writes all genes into one command file by default and prints
the number of genes and command files. It also returns:

```julia
out.genes
out.fitfile
out.commandfiles
```

## Submit On Biowulf

Submit the single command file with Biowulf `swarm`. For large gene panels, use
`-b` at submission time to bundle several gene commands into each Slurm array
task:

```bash
cd run-RamosNELFA_NEG_IFNa_rep1
swarm -f fit.swarm -b 20 -g 4 -t 2 --time 2:00:00 --merge-output --module julia
```

Here `-t 2` matches `nchains=2, nthreads=1`. The `-b 20` option makes each
Slurm array task run 20 gene commands serially; it reduces array size without
changing the resources used by each Julia command. Set `batchsize=<N>` in
`makeswarm_genes` only if you deliberately want several numbered command files.
Increase memory or walltime if your model settings require it.

## Local Or Other Schedulers

The generated `.swarm` file is a plain command list, so it can be run
sequentially from a regular shell:

```bash
bash fit.swarm
```

For Slurm or GNU Parallel wrappers, request a scheduler launcher:

```julia
makeswarm_genes(
    datapath = "data/HCT116_testdata",
    datacond = "MOCK",
    resultfolder = "HCT116_test",
    filedir = "run-HCT116-testdata-rna",
    scheduler = :slurm,
    scheduler_jobs = 100,
    slurm_mem = "8G",
    slurm_time = "02:00:00",
)
```

or:

```julia
makeswarm_genes(
    datapath = "data/HCT116_testdata",
    datacond = "MOCK",
    resultfolder = "HCT116_test",
    filedir = "run-HCT116-testdata-rna",
    scheduler = :parallel,
    scheduler_jobs = 16,
)
```

## Result Tables

After fits finish, assemble the legacy RNA result files into summary CSVs:

```julia
write_dataframes_only(
    "results/HCT116_test",
    "data/HCT116_testdata";
    datatype = "rna",
    rna_truncation = :legacy,
)
```

This calls `make_dataframes`, assembles raw `rates_*.txt`, `measures_*.txt`,
and `param-stats_*.txt` files when needed, and writes summary CSV files.
Existing fit-time `burst_*.txt` and `optimized_*.txt` files are assembled into
`burst_*.csv` and `optimized_*.csv` and joined into the summaries. Use the same
`rna_truncation` setting as the fit when observed RNA moments are recomputed.

For key-based folders containing `rates_<key>.txt`, `param-stats_<key>.txt`,
and `info_<key>.jld2`, use:

```julia
write_dataframes_only_key("results/HCT116_test"; datatype = "rna")
```

## Model Variants

The same input files can be used for different model structures:

```julia
for (G, R, transitions) in [
    (2, 0, ([1, 2], [2, 1])),
    (3, 0, ([1, 2], [2, 1], [2, 3], [3, 2])),
    (2, 1, ([1, 2], [2, 1])),
]
    fit(
        G = G,
        R = R,
        S = 0,
        insertstep = 1,
        transitions = transitions,
        datatype = "rna",
        datapath = "data/HCT116_testdata",
        gene = "MYC",
        cell = "HCT116",
        datacond = "MOCK",
        resultfolder = "HCT116_test",
    )
end
```

Rates are reported in inverse minutes. Halflife metadata, when used, are read
from the project `data/halflives` folder.
