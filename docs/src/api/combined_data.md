# v1.10 CombinedData API

StochasticGene v1.10 introduces a new `CombinedData` path for multimodal fits. The goal is to let each elementary modality keep its own loader and likelihood while `fit` combines the scalar log-likelihoods and WAIC pointwise predictions in a stable order.

This API is intended to replace new uses of legacy combined datatypes such as `"tracerna"` and `"rnadwelltime"` over time. Those legacy datatypes are still supported for compatibility, but new multimodal workflows should prefer tuple or vector `datatype` values.

## Single vs Combined Data

Single-modality fits can still use a string or symbol:

```julia
fit(;
    datatype = "rna",      # also accepts :rna
    datapath = "HBEC_smFISH",
    gene = "CANX",
    cell = "HBEC",
)
```

Combined fits use a tuple or vector of modality names:

```julia
fit(;
    datatype = (:rna, :dwelltime),
    datapath = (
        rna = "HBEC_smFISH",
        dwelltime = [
            "dwelltime/CANX_ON.csv",
            "dwelltime/CANX_OFF.csv",
        ],
    ),
    gene = "CANX",
    cell = "HBEC",
    dwell_specs = [
        (
            unit = 1,
            onstates = [Int[], Int[], [2, 3]],
            dttype = ["ON", "OFF"],
        ),
    ],
)
```

The modality order is canonicalized. For example, `(:dwelltime, :rna)` and `(:rna, :dwelltime)` both construct the same `CombinedData` key order. This keeps dispatch, likelihood assembly, output naming, and WAIC prediction order stable.

## Supported Modalities

The combined API currently recognizes these modality symbols:

- `:rna`
- `:trace`
- `:dwelltime`
- `:grid` is reserved for future work

The v1.10 likelihood stack has focused support for the combinations already used in the refactor:

- `(:rna, :trace)` for the split equivalent of legacy `"tracerna"`
- `(:rna, :dwelltime)` for the split equivalent of legacy `"rnadwelltime"`
- single-element combined tuples such as `(:rna,)` are useful internally and in tests

Legacy strings such as `"tracerna"`, `"rnadwelltime"`, and `"rnaonoff"` continue to load their legacy data structures.

## `datapath` Forms

For combined data, prefer a `NamedTuple` keyed by modality:

```julia
datapath = (
    rna = "HBEC_smFISH",
    dwelltime = [
        "dwelltime/CANX_ON.csv",
        "dwelltime/CANX_OFF.csv",
        "dwelltime/CANX_ONG.csv",
    ],
)
```

Each entry is resolved recursively under `root/data` when needed, so the example above resolves to paths under `joinpath(root, "data", ...)`.

For transition workflows, v1.10 also accepts the legacy positional layout for the two common combinations:

```julia
# Equivalent to the NamedTuple form for (:rna, :dwelltime)
datapath = (
    "HBEC_smFISH",
    "dwelltime/CANX_ON.csv",
    "dwelltime/CANX_OFF.csv",
)
```

The keyed form is clearer and is recommended for scripts, saved run specs, and examples.

## Output Names

`CombinedData` derives output labels and gene names from its elementary legs. If all legs share the same `label` and `gene`, output stems look like ordinary single-data fits:

```text
rates_FISH_CANX_3331_2.txt
measures_FISH_CANX_3331_2.txt
param-stats_FISH_CANX_3331_2.txt
```

If legs have different labels or genes, the stem joins the distinct values with `+`.

## Retired Legacy Arguments

The public keyword surface for new code is:

- `root`: project root
- `datapath`: data folder/file path, or a `NamedTuple` for combined data
- `label`: output/data label stem; if empty, `fit` builds one from datatype, cell, and condition
- `resultfolder`: output folder, usually resolved under `root/results`
- `trace_specs`: trace observation metadata
- `dwell_specs`: dwell-time observation metadata

The older `infolder` and `inlabel` style arguments are retired. `fit(; key=...)` still ignores an old `infolder` key if it appears in an existing run-spec file so older saved jobs do not immediately break, but new specs should not write it. Likewise, legacy `traceinfo` and `dttype` entries can be consumed for migration, then dropped from the persisted run spec; use `trace_specs` and `dwell_specs` going forward.

## Likelihood Behavior

For `CombinedData`, likelihood evaluation dispatches by modality. Each leg computes its own likelihood and pointwise log-prediction vector. The combined likelihood is the sum of leg likelihoods, and the WAIC vector is the concatenation of leg vectors in canonical modality order.

The same `Options` object created by `fit` is threaded into the combined likelihood stack. For MH this is `MHOptions`; for gradient-based inference this is `NUTSOptions` or `ADVIOptions`, including the selected `likelihood_executor` and any gradient checkpoint settings.
