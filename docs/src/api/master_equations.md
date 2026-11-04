# Master-Equation Solvers

These optional numerical methods operate on RNA probability distributions.
They do not replace the default steady-state RNA fitter, HMM transition solver,
or MCMC sampler.

All rates are per minute; times are minutes. Joint distributions use
`P[configuration, m + 1]` for RNA count `m`, so Julia's `vec(P)` has the
configuration index varying fastest. Finite generators have **column** sums
zero.

## Transient RNA Distributions

```julia
using StochasticGene

transitions = ([1, 2], [2, 1])
components = StochasticGene.MComponents(transitions, 2, 0, 60, 1.0, "")
rates = [0.3, 0.5, 8.0]  # k12, k21, production; decay is in components
problem = transient_master_problem(components, rates)
P0 = transient_master_initial([0.625, 0.375], 60)

P = transient_master_closure(problem, P0, 1.0; ode_steps=200)
histogram = transient_master_marginal(P)
```

`transient_master_problem` supports production that leaves the finite
configuration unchanged. For GR/GRSM ejection events, which create RNA **and**
change configuration, use `transient_master_problem_general(components, rates)`.
The original adapter continues to reject that case explicitly.
Both component adapters require single-transcript ejection and scalar linear
decay. They reject count-burst operators instead of treating multiple-transcript
events as single births. Direct construction can be used with an output window
of one bin; the component adapter needs at least two bins to recover decay.

For direct construction, `TransientMasterProblemAB(A, B, decay)` describes
nonproductive transitions `A` and single-transcript productive transitions `B`.
`A + B` must be a finite generator. The count equations are

```text
dP_m/dt = A*P_m + B*P_(m-1) + decay*((m+1)*P_(m+1) - m*P_m).
```

| Method | Production | Step control | Notes |
|:--|:--|:--|:--|
| `transient_master_strang` | Diagonal B | positional `n_steps` | Existing second-order split |
| `transient_master_strang_richardson` | Diagonal B | positional base steps | Extrapolated fourth order when window loss is negligible; may give negative entries |
| `transient_master_strang_purebd` | General B | positional `n_steps` | Death split plus RK4 production/configuration step |
| `transient_master_closure` | General B | `ode_steps` | RK4 multiplier ODE, no operator-splitting error |
| `transient_master_closure_taylor` | General B | `n_steps`, `R` | Explicit order-R multiplier integration |
| `transient_master_closure_exp_taylor` | General B | `n_steps`, `R` | Exponential treatment of the constant generator; Taylor treatment of the remaining term |

The closure methods use a generating-function composition/multiplier
representation. Their `M_out` is an **output window**, not a reflecting cap:
probability outside it is not added to the last bin or renormalized. `P0` can
have a larger count support than the output window. Any initial mass omitted
from `P0` is, of course, unavailable to the solver.

Splitting methods instead discard mass outside their represented count window
at intermediate steps. Death can otherwise return that mass to lower counts.
Check both count-window and time-step convergence; these are different errors.
Comparisons to the existing capped full CME require negligible cap probability
or a deliberately larger reference window.

All listed time integrators are fixed-step. Defaults are convenience choices,
not error-tolerance guarantees. Double step counts and compare results; inspect
nonnegativity and mass without silently clipping or normalizing away errors.
The new solvers reject nonfinite results, materially negative probabilities,
or a material increase over the supplied initial mass, asking for smaller
steps. Passing those checks does not establish accuracy. Richardson
extrapolation is deliberately not clipped and can return negative entries.
An exponential treatment of one term does not guarantee stability of the full
scheme at arbitrarily large steps. Taylor orders are restricted to 1 through
20. The new closure/state-decay methods use Float64 work arrays; they are not
gradient-based inference interfaces. No universal performance ranking is claimed.

## State-Dependent Decay

`TransientMasterProblemStateMu(A, B, decay_vector)` allows one decay rate per
finite configuration, including zero.

- `transient_master_strang_statewise(problem, P0, t, n_steps)` handles diagonal
  productive `B`. It rejects off-diagonal productive transitions rather than
  dropping them.
- `transient_master_strang_shifted(problem, P0, t, n_steps; mu0=minimum(decay_vector))`
  also handles general `B`. It splits off a reference scalar decay and a
  residual state-dependent death operator. Require
  `0 <= mu0 <= minimum(decay_vector)` so the residual is a valid death process.
  It uses a dense joint matrix exponential and is intended for moderate joint
  dimensions, not large coupled-state spaces.

## Optional Steady-State Block Solver

```julia
M = StochasticGene.make_mat_M(components, rates)
p = StochasticGene.steady_state_vector(M;
    solver=:closure_thomas, nT=components.nT)

h = StochasticGene.steady_state(M, components.nT, 2;
    steady_state_solver=:closure_thomas) # two alleles, existing convolution
```

The block-Thomas method exploits the count-level block-tridiagonal structure
of the **existing capped** generator. It preserves that boundary convention;
it is not the transient closure's uncapped output window. It validates matrix
layout and rejects unsupported structures, including nonstandard count jumps
and state-dependent decay. Singular elimination falls back to the existing
robust nullspace solver. The default solver is unchanged.

ForwardDiff through supported nonsingular elimination is tested. No reverse-mode
rule is added: retain `solver=:augmented` for the existing reverse-mode path.
Benchmark on the actual model before choosing a solver; dense finite-state
blocks can be expensive for large coupled configurations.

## API

```@docs
StochasticGene.TransientMasterProblem
TransientMasterProblemAB
TransientMasterProblemStateMu
transient_master_problem
transient_master_problem_general
transient_master_initial
transient_master_strang
transient_master_strang_richardson
transient_master_strang_purebd
transient_master_closure
transient_master_closure_taylor
transient_master_closure_exp_taylor
transient_master_strang_statewise
transient_master_strang_shifted
transient_master_marginal
StochasticGene.steady_state_vector
StochasticGene.normalized_nullspace
StochasticGene.normalized_nullspace_augmented
StochasticGene.normalized_nullspace_svd
StochasticGene.normalized_nullspace_qr
StochasticGene.pullback_normalized_nullspace_augmented
StochasticGene._augmented_nullspace_try
StochasticGene.normalized_nullspace_closure_thomas
StochasticGene.make_mat_M
```
