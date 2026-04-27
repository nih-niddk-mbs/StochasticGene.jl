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

To exercise the new inference/gradient stack, use the matrix runner:

```julia
include("test/fullstack/run_fits.jl")

run_fullstack_inference_matrix!()
```

The default matrix is meant to run for a few minutes and show progress. It uses
`samplesteps = 300`, `warmupsteps = 150`, `maxiter = 150`, `n_mc = 4`, and
`maxtime = 300.0`. The default expected-recovery pairs are:

- `rna:nuts_forwarddiff`
- `rna:advi_zygote`
- `rna:advi_forwarddiff`
- `trace:nuts_forwarddiff`
- `trace:nuts_finite`

Override any of those defaults when you want a longer or narrower check:

```julia
run_fullstack_inference_matrix!(
    cases="trace",
    specs="nuts_zygote,advi_zygote_trace",
    samplesteps=10,
    warmupsteps=5,
    maxiter=10,
)
```

When `cases` or `specs` is supplied, the runner uses the Cartesian product of
selected cases and specs. This is useful for experimental checks, but not every
combination is expected to recover rates.

Named inference specs are:

- `mh_fast`: MH with the fast likelihood stack.
- `nuts_forwarddiff`: NUTS with `gradient = :ForwardDiff`.
- `nuts_finite`: NUTS with central finite-difference gradients.
- `nuts_zygote`: NUTS with `gradient = :Zygote`.
- `advi_zygote`: ADVI with `gradient = :Zygote`.
- `advi_forwarddiff`: ADVI with `gradient = :ForwardDiff`.
- `advi_finite`: ADVI with `gradient = :finite` (implemented with ForwardDiff on the ELBO).
- `advi_zygote_trace`: ADVI with `zygote_trace = true` and checkpointing, intended only for short trace smoke checks.

For trace data, ADVI is currently an experimental plumbing check rather than an
expected recovery check. Plain `advi_zygote` intentionally falls back to the
ForwardDiff/finite ELBO gradient path unless `zygote_trace = true`; reverse-mode
Zygote through the trace HMM is memory-heavy and should be treated as a short-trace
diagnostic only.

From the shell, set `SG_FULLSTACK_INFERENCE_MATRIX=1` and optionally restrict cases/specs:

```bash
SG_FULLSTACK_INFERENCE_MATRIX=1 \
SG_FULLSTACK_CASES=rna,trace \
SG_FULLSTACK_INFERENCE_SPECS=nuts_forwarddiff,advi_zygote \
julia --project=. test/fullstack/run_fits.jl
```

Fit outputs are written to `test/fullstack/results/`, which is ignored by git.

The fixture truth is defined in `generate_data.jl`:

- `GENE1` RNA-like unit rates: `[0.30, 0.15, 0.60, 1.0, 0.10]`.
- `GENE2` RNA-like unit rates: `[0.25, 0.18, 0.55, 1.0, 0.12]`.
- Coupled tracejoint rates: the two unit rate blocks above plus coupling `-0.20`.

The fit drivers fit only selected parameters for speed:

- `rna`, `rnacount`, `trace`, and `dwelltime` fit rates 1 and 2.
- `tracejoint` fits only the coupling parameter.

After each fit, the driver prints a rate-recovery summary comparing median fitted
parameters to the corresponding fixture truth. This is a coarse smoke check, not
a statistical pass/fail test.
