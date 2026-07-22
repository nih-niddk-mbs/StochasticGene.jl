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
)
```

If `decayrate=-1.0` (the default), `fit` tries to look up gene-specific decay
information from `data/halflives` under `root`. If allele metadata exist under
`data/alleles`, `fit` can also use those to set `nalleles`.

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
    batchsize = 4800,
    maxtime = 3600.0,
    samplesteps = 1_000_000,
    propcv = 0.01,
    project = "/home/carsonc/github/StochasticGene.jl/",
)
```

`makeswarm_genes` prints the number of genes, the batch size, and the number of
batch command files. It also returns:

```julia
out.genes
out.fitfile
out.commandfiles
```

By default, folder scanning uses filenames only. To restrict to genes that pass
the older metadata gate, use:

```julia
filter_metadata = true
```

## Submit On Biowulf

If `batchsize` splits the run into numbered swarm files, submit each file:

```bash
cd run-RamosNELFA_NEG_IFNa_rep1
for f in fit_*.swarm; do
    swarm -f "$f" -g 4 -t 2 --time 2:00:00 --merge-output --module julia
done
```

Here `-t 2` matches `nchains=2, nthreads=1`. Increase memory or walltime if
your model settings require it.

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
)
```

This calls `make_dataframes`, assembles raw `rates_*.txt`, `measures_*.txt`,
and `param-stats_*.txt` files when needed, and writes summary CSV files.

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
