# Documentation Audit

This page records the current documentation/docstring audit for the v1.11 beta
line. It is aimed at maintainers and should be updated whenever public APIs or
workflow defaults change.

## Current Source Of Truth

- User workflows: `README.md`, `docs/src/getting_started.md`, and
  `docs/src/cluster_batch_workflows.md`.
- Public API behavior: docstrings in `src/*.jl`, especially `fit.jl`,
  `biowulf.jl`, `analysis.jl`, `io.jl`, and `stage.jl`.
- Exported symbols: the `export` block in `src/StochasticGene.jl`. The old
  unused `src/exports.jl` file has been retired so there is only one canonical
  export list.
- Generated website navigation: `docs/make.jl`.
- Do not edit generated files under `docs/build`; update `docs/src` and rebuild.

## Audited Workflows

### RNA histogram and scRNA sweeps

Status: current.

- `makeswarm_genes(; datapath, datacond, ...)` is the replacement for the older
  v0.7-style folder sweep.
- Folder scanning uses `checkgenes` by default
  (`filter_metadata=true`), preserving the halflife/allele metadata gate.
- `filter_metadata=false` scans filenames only.
- `batchsize=nothing` writes all gene commands into one command file by default.
- For large Biowulf runs, submit the single file with `swarm -b <bundle>` rather
  than splitting into many numbered command files.
- `batchsize=<N>` is still available only when a user deliberately wants several
  command files.

Representative submit command:

```bash
swarm -f fit.swarm -b 20 -g 24 -t 4 --time 25:00 --merge-output --module julia
```

### RNA result dataframe stack

Status: current.

- Legacy/RNA-style folders:
  - `make_dataframes(resultfolder, datapath; ...)`
  - `write_dataframes_only(resultfolder, datapath; ...)`
  - `write_dataframes(resultfolder, datapath; measure=:AIC, ...)`
- Key-based folders:
  - `make_dataframes_key(resultfolder; datapath=nothing, ...)`
  - `write_dataframes_only_key(resultfolder; datapath=nothing, ...)`

Docs now distinguish these paths explicitly.

### Key-based trace and shared analyses

Status: current for entry points; deeper examples can be expanded.

- `write_traces_key(folder)` should be used for key-based trace prediction.
- `write_ONOFFhistograms_key(folder)` is the key-aware ON/OFF histogram path.
- `write_correlation_functions_key(folder; ...)` is the modern key-aware
  theoretical correlation path and can also write generalized observable-pair
  correlations.
- `write_shared_group_measures_from_samples(resultfolder, key; ...)` is the
  post-fit shared/group diagnostic path when posterior samples are available.

### Retired or compatibility-only inputs

Status: current.

- `infolder` and `inlabel` are retired public inputs. Use `root`, `datapath`,
  `label`, and `resultfolder`.
- `fit(; key=...)` still tolerates old run-spec fields during migration.
- `traceinfo` and `dttype` remain accepted in legacy scripts, but new trace and
  dwell workflows should prefer `trace_specs` and `dwell_specs`.

## Known Documentation Risks

- The package exports many developer/test/benchmark helpers from
  `src/StochasticGene.jl`. Not all of those are intended as user-facing APIs.
  Avoid adding manual examples for test helpers unless they are deliberately
  promoted.
- `docs/src/api/index.md` is broad and partly hand-written. When new public
  APIs are added, prefer adding specific pages or concrete workflow docs rather
  than expanding that file with generic placeholder APIs.
- Analysis APIs should be documented only when the functions exist in `src`.
  Placeholder names such as `profile_likelihood`, `posterior_mean`, or
  `trace_plots` should not appear in the manual unless implemented.
- Shared-fit documentation should keep the rate/parameter distinction clear:
  rate files are complete likelihood-ready rate vectors; fitted parameters are
  the reduced parameter vector used by the sampler.

## Audit Checklist For Future Changes

1. Search docs and source docstrings for stale keywords:

   ```bash
   rg -n "infolder|inlabel|batchsize|filter_metadata|traceinfo|dttype|fit_rna" README.md docs/src src
   ```

2. Search for old Biowulf split-file guidance:

   ```bash
   rg -n "fit_\\*\\.swarm|for f in fit_|batch command files|batchsize = (1000|4800)" README.md docs/src
   ```

3. Verify Julia source parses after docstring edits:

   ```bash
   julia --startup-file=no -e 'Meta.parseall(read("src/biowulf.jl", String)); Meta.parseall(read("src/io.jl", String)); Meta.parseall(read("src/analysis.jl", String)); println("parse ok")'
   ```

4. If dependencies are available, build docs:

   ```bash
   julia --project=docs docs/make.jl
   ```

5. Never patch `docs/build` directly.
