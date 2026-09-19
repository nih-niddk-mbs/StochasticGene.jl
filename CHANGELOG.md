# Changelog

## 2.0.6

Correctness patch for the theoretical trace-centering postprocessor:

- `write_correlation_functions_centered` raised a `MethodError` instead of the
  documented `ArgumentError` when a correlation file did not cover the full
  requested trace window, because its diagnostic `@warn` used an invalid
  backslash line continuation that Julia parsed as a division operator.
- `write_correlation_functions_centered_folder` (and the folder form of
  `write_correlation_functions_centered`) silently skipped existing
  `crosscorrelation-global_*` files instead of trace-centering them.
- Fixed a `MethodError` in the legacy trace-centering path
  (`_trace_center_legacy_correlation!`) caused by building its column/mean-product
  lookup as a tuple of pairs instead of a `Dict`.
- Unified the column convention between the legacy and generalized correlation
  CSVs: `cc_ON`/`cc`, `ac1_ON`/`ac_x`, etc. now always hold the raw,
  uncentered moment (unchanged meaning in every file), and a `cc_ON_centered`/
  `cc_centered`, `ac1_ON_centered`/`ac_x_centered`, etc. companion column is
  always present holding the centered covariance — a simple mean-product
  subtraction in base/global-centered files, or the full finite-window
  trace-centering result in trace-centered files. Previously the legacy schema
  conflated the two, storing a centered-plus-mean-product hybrid directly in
  `cc_ON` for trace-centered output.

## 2.0.5

Correctness patch for legacy hierarchical trace models:

- Evaluate individual hierarchy distributions entirely in physical parameter
  space, after the inverse sampling transforms have been applied.
- Interpret the first hierarchical block as population means and the second as
  positive coefficients of variation. Positive individual parameters use a
  log-normal conditional distribution parameterized by arithmetic mean and CV;
  signed individual parameters use a normal conditional distribution.
- Include the individual inverse-transform Jacobian exactly once in the
  posterior density sampled by MH, NUTS, and ADVI.
- Keep HMM methods responsible only for the observation likelihood. Hierarchy
  terms are now assembled by the higher-level trace likelihood, including grid
  models, while pointwise trace terms remain data-only for WAIC.

This corrects mixed-coordinate hierarchy calculations present from version
1.5.0 through 2.0.4. Hierarchical posterior values and fitted results can change.
Existing rate and info files are not rewritten; continuing an old run uses its
saved physical rates as an initial condition for the corrected posterior.

## 2.0.4

ON/OFF dwell-time output improvement:

- Write the directly calculated `ON_CDF` and `OFF_CDF` alongside the existing
  normalized PDF columns.
- Add uniform `bin_size` and `maxtime` options, in minutes, to the key-based and
  legacy folder generators. Use `bin_size=100/60` for 100-second frames.
- Include `SimON_CDF` and `SimOFF_CDF` when simulation output is requested.

Existing `ON` and `OFF` PDF columns and direct explicit-bin calls are preserved.

## 2.0.3

Correctness and documentation patch for coupled simulation:

- Make coupled simulation use the same connection semantics as fitting for
  multiple target transitions, reciprocal coupling, hidden three-unit models,
  and repeated entries in the unit-model map.
- Support both RNA reporter coupling definitions: `Rany` applies one coupling
  factor when any reporter position is occupied, whereas `Rsum` adds one
  contribution for every occupied reporter position.
- Generalize stationary-state indexing beyond two units and preserve coupled
  `a_grid` settings in trace simulation.
- Validate coupled model maps, rate/noise vector lengths, and coupling sign
  metadata instead of silently interpreting malformed inputs.
- Replace obsolete coupled-simulator examples and document the canonical
  `(unit_model, connections[, sign_modes])` representation and rate layout.

These changes can alter simulated trajectories from coupled models that use
multiple connections, `Rany`/`Rsum`, reciprocal coupling, hidden units, or a
reused model definition. Fit output formats are unchanged.

## 2.0.2

Correctness patch for inhibitory and free `Rsum` coupled models:

- Constrain a tied `Rsum` coupling shared by `m` simultaneously occupiable R
  positions to `γ > -1/m`, rather than constraining each expanded connection
  independently to `γ > -1`.
- Apply the overlap-dependent bound in both ordinary and hierarchical fits,
  including fixed-effect detection of tied R-position connections.
- Keep default inhibitory `Rsum` priors strictly inside the valid interval for
  models with any number of R steps.
- Preserve the ordinary `γ > -1` bound for `Rany` and for tied couplings whose
  source states are mutually exclusive.

Existing `Rsum` rate files fitted with earlier releases are not rewritten
automatically. If a historical fitted coupling has `1 + mγ <= 0`, truncate or
otherwise migrate the starting value into the open interval before continuing
the fit. Explicit `Rsumk` and `Ranyk` CSV tokens can be used to distinguish the
two coupling definitions in new workflows.

## 2.0.1

Correctness patch for coupled simulation and correlation analysis:

- Correct theoretical correlation centering, including unit-2
  autocorrelation centering, so centered correlations decay to zero.
- Correct coupled simulation when multiple connections affect the same target
  transition.
- Correct event-count and physical-time burn-in sequencing so recorded traces
  begin from the equilibrium-time distribution rather than a reaction epoch.
- Correct sparse simulated trace endpoints.
- Add `nexperiments` to `simulate_trials` for threaded, repeated-experiment
  ON-ON correlation validation and finite-experiment intervals.

## 2.0.0

StochasticGene 2.0 consolidates the inference, shared-parameter, analysis, and
batch-generation work developed on the 1.11 release line. See the
[2.0 release notes](docs/src/release_notes.md) for highlights and migration
guidance.

## 1.11.0

Beta release for the 2.0 inference and workflow stack.

## 1.10.1

Previous stable release.
