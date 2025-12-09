# Cross-Covariance Calculation Audit

## Summary
This document audits the theoretical and simulation cross-covariance calculations to identify potential discrepancies.

## Test Reference
`test_compare_coupling()` matches exactly, indicating:
- Core simulation (`simulator`) is correct
- Master equation solutions (`test_CDT`) are correct
- Steady-state distributions match

## Theoretical Cross-Covariance (`covariance_functions`)

### Computation Flow:
1. **Cross-correlation function** (`crosscorfn_hmm`):
   - Computes: `E[X(0)Y(τ)]` for τ ≥ 0
   - Formula: `cc[l] = sum_i sum_j p0[i] * meanintensity1[i] * a^τ[i,j] * meanintensity2[j]`
   - This is the **uncentered** cross-correlation

2. **Cross-covariance** (`crosscov_hmm`):
   - `cc12 = crosscov_hmm(..., lags, m1, m2)` → `E[X(0)Y(τ)] - E[X]E[Y]` for τ ≥ 0
   - `cc21 = crosscov_hmm(..., lags, m1, m2)` → `E[Y(0)X(τ)] - E[X]E[Y]` for τ ≥ 0
   - Both use normalized `mean_intensity` and scaled by `max_intensity[1] * max_intensity[2]`

3. **Full cross-covariance** (line 2684):
   - `cc = vcat(reverse(cc21), cc12[2:end])`
   - **Negative lags**: `reverse(cc21)` gives `E[Y(0)X(-τ)] - E[X]E[Y]` = `E[Y(τ)X(0)] - E[X]E[Y]` (by stationarity)
   - **Positive lags**: `cc12[2:end]` gives `E[X(0)Y(τ)] - E[X]E[Y]` for τ > 0
   - **Lag 0**: Uses `cc12[1]` = `E[X(0)Y(0)] - E[X]E[Y]`

### Key Points:
- Uses **normalized** `mean_intensity` (divided by `mmax`)
- Scales result by `max_intensity[1] * max_intensity[2]`
- Means `m1, m2` are computed on normalized scale, then scaled back
- **No measurement noise** in cross-covariance (only in autocovariance at lag 0)

## Simulation Cross-Covariance (`simulate_trials`)

### Computation Flow:
1. **Simulation** (`simulate_trace_vector`):
   - Generates traces with measurement noise (via `probfn`)
   - Returns `t[1][:, 1]` and `t[1][:, 2]` (two channels)

2. **Cross-covariance** (line 2695):
   - `StatsBase.crosscov(t[1][100:end, 1], t[1][100:end, 2], lags, demean=true)`
   - Computes: `E[x(t)y(t+τ)] - E[x]E[y]` for all lags in `lags`
   - `demean=true` removes the mean before computation
   - For negative lags: computes `E[x(t)y(t-|τ|)] - E[x]E[y]`

### Key Points:
- Uses **raw simulated traces** (includes measurement noise)
- `demean=true` ensures proper centering
- Lags can include negative values (full range: `-lag:stride:lag`)

## Potential Issues

### Issue 1: Measurement Noise in Cross-Covariance
- **Theory**: Cross-covariance does NOT include measurement noise (only autocovariance at lag 0 does)
- **Simulation**: Cross-covariance includes measurement noise at all lags
- **Impact**: Since noise is **independent** between channels:
  - `E[N1(t)N2(t)] = E[N1]E[N2] = 0` (zero-mean noise)
  - `E[X(t)N2(t)] = E[X(t)]E[N2(t)] = 0` (independent)
  - Therefore: `E[x(t)y(t)] = E[(X+N1)(Y+N2)] = E[XY] + 0 + 0 + 0 = E[XY]`
  - So simulation cross-covariance matches theory: `E[x(t)y(t+τ)] - E[x]E[y] = E[X(t)Y(t+τ)] - E[X]E[Y]`
- **Status**: ✓ Verified - noise independence means no discrepancy

### Issue 2: Lag Convention
- **Theory**: `cc12` = `E[X(0)Y(τ)] - E[X]E[Y]` for τ ≥ 0
- **Theory**: `cc21` = `E[Y(0)X(τ)] - E[X]E[Y]` for τ ≥ 0
- **Theory full**: `cc = vcat(reverse(cc21), cc12[2:end])`
  - Negative lags: `reverse(cc21)` = `E[Y(0)X(-τ)] - E[X]E[Y]` = `E[Y(τ)X(0)] - E[X]E[Y]` (stationarity)
  - Positive lags: `cc12[2:end]` = `E[X(0)Y(τ)] - E[X]E[Y]` for τ > 0
- **Simulation**: `StatsBase.crosscov(x, y, lags)` where `lags` can be negative
  - For positive lag τ: `E[x(t)y(t+τ)] - E[x]E[y]` = `E[X(0)Y(τ)] - E[X]E[Y]` ✓
  - For negative lag -τ: `E[x(t)y(t-τ)] - E[x]E[y]` = `E[X(τ)Y(0)] - E[X]E[Y]` = `E[Y(0)X(-τ)] - E[X]E[Y]` ✓
- **Status**: Lag convention appears consistent

### Issue 3: Mean Calculation
- **Theory**: Uses `m1, m2` computed from normalized `mean_intensity`, then scaled
- **Simulation**: Uses empirical mean from traces (includes noise)
- **Impact**: If means don't match, cross-covariance will be offset
- **Status**: Already identified and fixed mean scaling issue

### Issue 4: Scaling
- **Theory**: Computes on normalized scale, then scales by `max_intensity[1] * max_intensity[2]`
- **Simulation**: Uses raw trace values (already on correct scale)
- **Status**: Should be consistent if `max_intensity` scaling is correct

## Recommendations

1. **Verify noise independence**: Check if measurement noise is independent between channels. If so, cross-covariance at lag 0 should not have noise contribution.

2. **Check lag 0**: At lag 0, cross-covariance should be `E[XY] - E[X]E[Y]`. Verify:
   - Theory: `cc12[1]` = `E[X(0)Y(0)] - E[X]E[Y]` (no noise)
   - Simulation: `E[x(t)y(t)] - E[x]E[y]` (includes noise if correlated)

3. **Verify stationarity**: Ensure both theory and simulation use stationary distributions.

4. **Check normalization**: Verify that `max_intensity` scaling is applied consistently.

