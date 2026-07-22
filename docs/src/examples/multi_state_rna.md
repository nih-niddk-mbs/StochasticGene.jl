# Multi-State RNA Analysis

This example shows how to fit `datatype="rna"` histograms with more than two
gene states.

## Data Layout

The RNA loader uses the same `GENE_COND.txt` histogram layout as the basic RNA
example:

```text
data/HCT116_testdata/
├── CENPL_MOCK.txt
└── MYC_MOCK.txt
```

The first column is the RNA copy-number histogram.

## Three-State Model

A common three-state topology has reversible transitions between states 1 and 2
and between states 2 and 3:

```julia
using StochasticGene

transitions = ([1, 2], [2, 1], [2, 3], [3, 2])

fits, stats, measures, data, model, options = fit(
    G = 3,
    R = 0,
    S = 0,
    insertstep = 1,
    transitions = transitions,
    datatype = "rna",
    datapath = "data/HCT116_testdata",
    gene = "MYC",
    cell = "HCT116",
    datacond = "MOCK",
    resultfolder = "HCT116_test",
    fittedparam = [1, 2, 3, 4],
    nchains = 2,
    maxtime = 60.0,
    samplesteps = 1_000_000,
    propcv = 0.01,
)
```

For RNA-only `G` models with `R=0`, the fitted parameters above correspond to
the transition rates. The final rate is the RNA decay rate; with
`decayrate=-1.0` (default), `fit` tries to use gene-specific halflife metadata
when available.

## Model Comparison

You can fit multiple state topologies to the same histogram and compare the
written `measures_*.txt` / summary CSV outputs:

```julia
models = [
    (G = 2, transitions = ([1, 2], [2, 1])),
    (G = 3, transitions = ([1, 2], [2, 1], [2, 3], [3, 2])),
]

for spec in models
    fit(
        G = spec.G,
        R = 0,
        S = 0,
        insertstep = 1,
        transitions = spec.transitions,
        datatype = "rna",
        datapath = "data/HCT116_testdata",
        gene = "MYC",
        cell = "HCT116",
        datacond = "MOCK",
        resultfolder = "HCT116_test",
        nchains = 2,
        maxtime = 60.0,
        samplesteps = 1_000_000,
    )
end
```

Then assemble result tables:

```julia
write_dataframes_only(
    "results/HCT116_test",
    "data/HCT116_testdata";
    datatype = "rna",
)
```

## Gene Panels

The same model can be run over every matching gene file with
`makeswarm_genes`:

```julia
out = makeswarm_genes(
    datapath = "data/HCT116_testdata",
    datacond = "MOCK",
    cell = "HCT116",
    resultfolder = "HCT116_test",
    filedir = "run-HCT116-three-state",
    G = 3,
    R = 0,
    S = 0,
    insertstep = 1,
    transitions = ([1, 2], [2, 1], [2, 3], [3, 2]),
    fittedparam = [1, 2, 3, 4],
    nchains = 2,
    nthreads = 1,
    batchsize = 1000,
)
```

`makeswarm_genes` scans filenames by default. Use `filter_metadata=true` only
when you want to restrict the gene list through `checkgenes`.
