# Analysis Functions

This page documents the current post-fit analysis surface. It intentionally
lists functions that exist in the package today; older placeholder names such as
`posterior_mean`, `trace_plots`, `profile_likelihood`, and `cross_validation`
are not public StochasticGene APIs.

Most workflows start from files written by [`fit`](fit.md), usually
`rates_*.txt`, `param-stats_*.txt`, `measures_*.txt`, and, for key-based runs,
`info_<key>.jld2`.

## Result Tables

### RNA-style result folders

Use `make_dataframes` when the result folder uses the older RNA-style naming
scheme:

```julia
dfs = make_dataframes(
    "results/HCT116_test",
    "data/HCT116_testdata";
    datatype = "rna",
)
```

Use `write_dataframes_only` to write the returned summaries, or
`write_dataframes` to also write winner tables:

```julia
write_dataframes_only(
    "results/HCT116_test",
    "data/HCT116_testdata";
    datatype = "rna",
)

write_dataframes(
    "results/HCT116_test",
    "data/HCT116_testdata";
    datatype = "rna",
    measure = :AIC,
)
```

These helpers assemble `rates_*.txt`, `param-stats_*.txt`, and
`measures_*.txt` files when `assemble=true`, then add observed RNA moments for
`datatype="rna"` or `datatype="rnacount"`.

### Key-based result folders

Use the key-specific helpers when files are named like `rates_<key>.txt` and
the folder contains `info_<key>.jld2`:

```julia
dfs = make_dataframes_key("results/HCT116_test"; datatype = "rna")
write_dataframes_only_key("results/HCT116_test"; datatype = "rna")
```

`make_dataframes_key` preserves the run key in a `Key` column and reads
`info_<key>.jld2` when present to add `Gene`, `Condition`, `Model`, `Nalleles`,
and `DataPath` columns.

## Prediction Helpers

`predictedarray` and `predictedfn` are lower-level likelihood helpers used by
the fitting stack. They generate model predictions for a rate vector, data
object, and model object:

```julia
fits, stats, measures, data, model, options = fit(...)
r = get_rates(stats.medparam, model, false)
pred = predictedarray(r, data, model)
```

Use these only after `fit` has constructed compatible `data` and `model`
objects. For ordinary users, the file-writing functions below are usually more
convenient.

## Trace Predictions

Use `write_traces_key` for modern key-based trace folders:

```julia
write_traces_key("results/5Prime-CRISPR-2026-07-22/")
```

This walks `info_*.jld2`, loads the matching `rates_<key>.txt`, and writes
predicted trace tables. It handles individual and shared-fit outputs by using
the run-spec metadata and the complete rate files associated with each group.

The lower-level `make_traces_dataframe` methods convert model-predicted traces
to `DataFrame`s after the caller has already supplied compatible rates, data,
and model metadata.

## ON/OFF And Residency Outputs

For key-based runs:

```julia
write_ONOFFhistograms_key("results/my-run")
write_joint_residence_prob_onoff_key("results/my-run")
write_residency_G_folder("results/my-run")
```

`write_ONOFFhistograms_key` is the key-aware replacement for older
`write_ONOFFhistograms` calls when the folder is organized around
`info_<key>.jld2`.

`write_joint_residence_prob_onoff_key` computes ON/OFF joint residence
probabilities from key metadata and matching rate files.

`write_residency_G_folder` computes G-state residence outputs for rate files in
a folder. For key-based folders, prefer keeping the `info_<key>.jld2` files with
the rates so the analysis code can recover model metadata without guessing.

## Correlation Functions

`write_correlation_functions_key` is the modern key-aware entry point for
theoretical correlation functions:

```julia
write_correlation_functions_key(
    "results/coupled-run";
    lags = collect(0:200),
)
```

By default it writes the historical ON/reporter cross-correlation CSV and the
generalized correlation CSV. The generalized path can evaluate arbitrary
observable pairs supported by `correlation_observable`, including ON indicators
and reporter counts.

For repeated correlation calculations, use:

```julia
ctx = build_correlation_context(r, transitions, G, R, S, insertstep, probfn, coupling, lags)
res = correlate_observables(ctx, :ON1, :ON2)
```

That avoids rebuilding the expensive model/correlation context for every
observable pair.

## Burst Prediction Analyses

Two AUC-style helpers address burst observability and identifiability:

```julia
known = simulate_burst_observability_auc(r, transitions, G, R, S, insertstep)

fitcheck = simulate_fit_burst_identifiability_auc(
    r,
    transitions,
    G,
    R,
    S,
    insertstep;
    fit_kwargs = Dict(:maxtime => 600.0, :nchains => 4),
)
```

`simulate_burst_observability_auc` asks: given fixed rates, how well can the HMM
recover simulated ground-truth bursts from noisy traces?

`simulate_fit_burst_identifiability_auc` asks: if traces are simulated from
known rates and then refit through the normal `fit` stack, how well do the
fitted rates recover the same burst labels?

## Shared-Fit Group Diagnostics

When posterior samples are saved for a shared fit, use
`write_shared_group_measures_from_samples` to recompute measures by group:

```julia
write_shared_group_measures_from_samples(
    "results/shared-CRISPR-2026-07-20",
    "shared-CRISPR-3301-gene",
)
```

This is intended for diagnosing whether a shared fit is favoring one data group
over another. It uses saved thinned posterior samples and the existing
likelihood machinery rather than a separate fitting path.

## See Also

- [Model fitting](fit.md)
- [Data loading](load_data.md)
- [Utility functions](utilities.md)
- [Cluster and batch workflows](../cluster_batch_workflows.md)
