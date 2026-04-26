# Developer Full-Stack Fixtures

This folder is for developer-only full-stack checks that read data files and write
results. It is intentionally not included from `test/runtests.jl`.

Generate the small fixture data with:

```bash
julia --project=. test/fullstack/generate_data.jl
```

Fixture convention:

- `genes = ("GENE1", "GENE2")`
- `cell = "CELL"`
- `datacond = "COND"`

The generator writes the current loader formats:

- `data/rna/GENE*_COND.txt`: one-column RNA histograms.
- `data/rnacount/GENE*_COND.txt`: two columns, RNA count and per-cell yield factor.
- `data/trace/*_GENE*_COND.txt`: paired two-unit trace files with `(frame, time, intensity)` columns; trace loading reads column 3.
- `data/dwelltime/GENE*_COND_OFF.txt`: two columns, bins and dwell-time histogram values.

Use the same trace folder for a single-trace fit and a joint trace fit:

- Single trace: `datatype = "trace"`, `gene = "GENE1"`, `datacond = "GENE1_COND"`.
- Joint trace: `datatype = "tracejoint"`, `datacond = ["GENE1_COND", "GENE2_COND"]`.

The scalar trace loader filters by the trace label (`datacond`), so the single-trace
case uses `GENE1_COND` to select only one gene from the shared trace folder.

Run the smoke-fit driver with:

```bash
julia --project=. test/fullstack/run_fits.jl
```

To run only selected cases:

```bash
SG_FULLSTACK_CASES=rna,trace julia --project=. test/fullstack/run_fits.jl
```

From Julia, you can include the driver and override sampling options:

```julia
include("test/fullstack/run_fits.jl")
run_fullstack_fits!(cases="rna,tracejoint", samplesteps=5, maxtime=30.0)
```

Fit outputs are written to `test/fullstack/results/`, which is ignored by git.
