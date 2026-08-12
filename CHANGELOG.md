# Changelog

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
