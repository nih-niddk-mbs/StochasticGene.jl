# Changelog

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
