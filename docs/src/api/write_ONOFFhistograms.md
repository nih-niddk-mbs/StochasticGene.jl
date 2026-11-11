# ON/OFF dwell-time outputs

`write_ONOFFhistograms` and `write_ONOFFhistograms_key` calculate theoretical
ON/OFF dwell-time distributions from fitted model rates. The key-based method
is recommended for current result folders because it reads model metadata from
`info_<key>.jld2` rather than inferring the model from a filename.

## Key-based use

```julia
using StochasticGene

write_ONOFFhistograms_key(
    "results/my-run";
    ratetype = "median",
    bin_size = 100 / 60,
    maxtime = 200.0,
)
```

StochasticGene rates and times use minutes, so a 100-second imaging interval is
`100 / 60` minutes. The call above evaluates the distribution at
`100/60, 200/60, ...` through 200 minutes.

To use nonuniform or precomputed time points, pass them directly:

```julia
write_ONOFFhistograms_key(
    "results/my-run";
    bins = collect(0.5:0.5:100.0),
)
```

An explicit `bins` vector takes precedence over `bin_size` and `maxtime`.

Bins are **positive right edges**, not centers. For a uniform grid
`dt:dt:T`, bin `j` is `((j-1)dt, j*dt]`. Zero-start grids are rejected;
check the original data's bin definitions rather than shifting labels blindly.
The low-level `dwelltimePDF` and `dwelltimeCDF` functions accept an explicit
`leftedge` for offset/nonuniform grids or a single bin. Otherwise the first
left edge is inferred from the first spacing, bounded below by zero.

## Output columns

Each `ONOFF_<key>.csv` contains:

| Column | Meaning |
|---|---|
| `time` | Dwell time in minutes |
| `ON` | Historical normalized ON PDF bin probability |
| `OFF` | Historical normalized OFF PDF bin probability |
| `ON_CDF` | Directly calculated ON cumulative probability |
| `OFF_CDF` | Directly calculated OFF cumulative probability |

The CDF is computed by the master-equation calculation and the PDF is derived
from its finite differences. The CDF columns are therefore preferable when a
plot or comparison requires cumulative probabilities; callers no longer need
to reconstruct them from the normalized PDF columns.

Despite their historical names, `ON` and `OFF` are **bin probability masses**,
normalized over the represented time window, not density heights. Their
cumulative sums need not equal the raw `ON_CDF`/`OFF_CDF` columns when the
window omits probability mass. This conditioning describes truncation; it does
not account for right-censored events or an omitted overflow bin. Dwell-time
likelihood input values are counts or deliberately selected likelihood weights.

With `simulate=true`, the output additionally contains `SimON`, `SimOFF`,
`SimON_CDF`, and `SimOFF_CDF`.

## Direct rate-vector use

```julia
bins = collect((100 / 60):(100 / 60):200.0)

df = write_ONOFFhistograms(
    rates,
    transitions,
    G,
    R,
    S,
    insertstep,
    bins;
    outfile = "ONOFF_model.csv",
)
```

The direct method requires a complete model rate vector in the standard rate
order. It returns the same `DataFrame` that it writes.

## Folder compatibility method

Legacy rate-file folders can still be processed with:

```julia
write_ONOFFhistograms(
    "results/legacy-run";
    bin_size = 100 / 60,
    maxtime = 200.0,
)
```

For key-based folders, prefer `write_ONOFFhistograms_key`.
