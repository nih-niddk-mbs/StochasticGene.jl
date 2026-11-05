# Choose a Model and Prepare Data

This page helps you decide which function to call, which model to start with, and
how your files must look. For definitions of the states and steps, see the
[Package overview](package_overview.md).

## What data do you have?

| You have | Call | Notes |
|----------|------|-------|
| Counts of RNA per cell (smFISH, scRNA-seq) as a histogram | [`fit_rna`](api/fit.md#Quick-start-helpers) or `fit(datatype = "rna")` | Steady-state data. Constrains the model less than traces. |
| Per-cell RNA counts with a yield factor per cell | `fit(datatype = "rnacount")` | Two columns: count and yield factor. |
| Live-cell intensity traces (MS2 or similar) | [`fit_traces`](api/fit.md#Quick-start-helpers) or `fit(datatype = "trace")` | Needs a reporter, so `R >= 1`. |
| ON/OFF dwell-time histograms, with or without RNA | `fit(datatype = (:rna, :dwelltime))` | Multimodal fit; see [CombinedData](api/combined_data.md). |
| Traces from two coupled signals | `fit(datatype = "tracejoint")` | See [Joint Trace Analysis](examples/joint_trace.md). |

Several data types can be fit together with a tuple such as `(:rna, :trace)`; each
data type contributes its own likelihood and they are summed.

## Which model should I start with?

The model is set by four numbers:

- `G`: number of gene states (2 gives the usual ON/OFF telegraph model).
- `R`: number of steps along the transcribed region that the polymerase passes.
  The last step releases the RNA.
- `S`: number of splice sites, at most `R - insertstep + 1` (larger values are reduced).
- `insertstep`: the step at which the reporter sequence appears (1 is the first step).

Starting points:

| Question | Start with |
|----------|------------|
| RNA histogram only | `G = 2`, `R = 0` |
| Live-cell traces | `G = 2`, `R = 1` or more, `S = 0` |
| Does a third gene state help? | Compare `G = 2` and `G = 3` on the same data |
| Splicing data | Add `S` only when the data can see it |

Add states and steps one at a time, and keep the extra ones only if they improve
the fit enough to matter. More states than the data can support give very wide
posteriors rather than better answers. Compare models with the WAIC and AIC that
`fit` prints.

For `G >= 3` the default transitions form a reversible chain
(`1 <-> 2 <-> 3`). Pass `transitions` to choose another layout such as a cycle.

## Units

- Rates are per **minute**.
- Frame intervals (`interval` in `fit_traces`, `trace_specs`) are in **minutes**.
- The `halflife` argument of `fit_rna` and `fit_traces` is in **minutes**.
- Half-lives in `data/halflives/<cell>_halflife.csv`, used by `fit` when no
  `decayrate` is given, are in **hours**. Do not mix these two up.

Without a half-life, the quick-start functions fix the decay rate at 1.0 per
minute. Other rates are then only meaningful relative to decay.

## File formats

### RNA histogram (`"rna"`)

- One text or CSV file per gene and condition, named `<gene>_<datacond>.txt`, in the
  `datapath` folder.
- The first column is the histogram: row `i` is the number of cells with `i - 1` RNA.
- `rna_histogram(counts)` builds this from per-cell counts.
- `rna_truncation = :legacy` (default) trims long tails; use `:none` to keep all bins.

### Per-cell counts (`"rnacount"`)

- Same file name. Column 1 is the RNA count of each cell, column 2 is that cell's
  yield factor.

### Traces (`"trace"`)

- Text or CSV files, one file per trace, one frame per row.
- Intensity is read from column `datacol` (default 3, as in `frame, time, intensity`).
- A file is used when its name contains the `datacond` string. Folders are searched
  recursively.
- Observation window and frame interval go in `trace_specs`
  (`unit`, `interval`, `start`, `t_end`, `zeromedian`, ...). `fit_traces` sets these for you.
- By default (`zeromedian = false`) `fit_traces` uses the raw intensities and estimates the
  `noisepriors` from the pooled data (a two-component Gaussian fit; reliable for `R = 1`, rough
  for larger `R`), and logs the values. Pass your own `noisepriors` if you know the noise.
- With `zeromedian = true` each trace is median-subtracted and scaled, and `noisepriors`
  default to `[0, 0.1, 1, 0.1]` on that normalized scale. If the real noise differs from the
  noise priors the fitted rates are biased: in a simulated test the OFF-to-ON rate came out about
  30 to 100 times too small. Prefer raw intensities.
- `noisepriors` are, in order, background mean, background SD, signal mean and
  signal SD.

### Dwell times

- Two columns: bin values and histogram counts, one file per ON or OFF type, for
  example `GENE_COND_OFF.txt`.

### Optional metadata

- `data/alleles/<cell>_alleles.csv`: allele counts per gene. If absent, `nalleles`
  (default 1) is used.
- `data/halflives/<cell>_halflife.csv`: half-lives in hours.
- `rna_setup()` creates the folder layout and copies example files.

## Run length and convergence

- RNA histogram fits are cheap. The default 100,000 samples is a reasonable start.
- Trace fits are slower and need a long adaptive warmup to mix well. `fit_traces` sets
  `warmupsteps` equal to `samplesteps` (default 40,000 each). In a simulated test, a fit with no
  warmup got stuck with the switching rates near zero and a burst size about 60 times too large.
  Use many more samples for real data.
- Check the diagnostics printed after sampling: R-hat close to 1 (below about 1.05)
  and a large effective sample size. A large R-hat means the chains disagree and the
  result should not be trusted. Run longer or more chains.
- Switching and initiation rates can be weakly identified from few traces, and a random-walk
  sampler can settle in a wrong region. In a simulated test (30 traces, 300 frames) the fit came
  back with a burst size of about 4.7 against a true 2.5 even with the noise correct, and the
  true rates had a higher likelihood than anything the sampler found. Use several chains from
  different starts, compare their likelihoods, and require R-hat near 1 and a large effective
  sample size before trusting a result. Report ranges, not single values.
- The posterior median is the recommended summary. See [the fit reference](api/fit.md#Priors,-indices,-outputs)
  for what the rows of `rates_*.txt` mean.
