# Why cross-correlated proposals (propcv from file) may not work well

This note traces the flow of **matrix/tuple propcv** (proposal covariance from a previous run) through `fit.jl` and `metropolis_hastings.jl` and lists likely reasons it "never quite worked very well."

---

## 1. Flow summary

- **fit.jl**: `propcv < 0` (scalar) → `get_propcv()` reads param-stats file → returns `(cov_matrix, abs(propcv))` or falls back to scalar.
- **make_structures** → **load_model**: tuple is stored in `model.proposal`.
- **run_mh**: `initial_proposal(model)` builds `d = proposal_dist(param, model.proposal, model)`.
- **proposal_dist(param, cv::Tuple, model)** calls `proposal_dist(param, cv[1], model, cv[2])` → **proposal_dist(param, cv::Matrix, model, indiv)** → `MvNormal(param, cv)` (or hierarchical split).
- Parameter space: `param = get_param(model)` is in **transformed** space (e.g. log). The file’s **last** block is `covlogparam` = covariance of **transformed** samples, so the read matrix is in the same space as `param`. ✓

---

## 2. Likely causes

### A. Warmup overwrites your covariance

**Where:** `metropolis_hastings.jl` → `run_mh` → `warmup`, then `sample(..., proposalcv, ...)`.

When `warmupsteps > 0` and adaptation runs (step > 1000 and accept rate > 0.15), the code **replaces** the initial proposal with the empirical covariance from the warmup run:

```julia
proposalcv = covparam * scaling   # or diag_cov fallback
d = proposal_dist(param, proposalcv, model)
return ..., proposalcv, ...
```

So a carefully chosen covariance loaded from file is **only used during warmup**; after that it is discarded and replaced by the warmup estimate. If warmup is short or has low/high acceptance, that estimate can be noisy, ill-conditioned, or poorly scaled, so the “cross-correlated” setup never actually drives the main sampling.

**Suggestion:** Add an option (e.g. “keep initial proposal” or “don’t adapt when proposal is a matrix/tuple”) so that when the user explicitly provides a covariance (e.g. from file), warmup does not overwrite it. Alternatively, only adapt when the initial proposal is scalar/vector, not matrix/tuple.

---

### B. No `proposal_dist` for `Diagonal` → possible error in warmup fallback

**Where:** `metropolis_hastings.jl` → `warmup`: when the empirical covariance is not positive definite, the code falls back to a diagonal matrix:

```julia
diag_cov = Diagonal(diag(covparam) * scaling)
proposalcv = diag_cov
d = proposal_dist(param, diag_cov, model)
```

There is **no** `proposal_dist(param, cv::Diagonal, model)` (only `Matrix` and scalar/vector/tuple). So this call can **method-error** when the warmup covariance is not PSD and the fallback is used.

**Suggestion:** Add a method that accepts a diagonal covariance, e.g.:

- `proposal_dist(param::Vector, cv::Diagonal, model, indiv=0.001) = proposal_dist(param, Matrix(cv), model, indiv)`, or
- implement it directly (e.g. build a product of Normals from the diagonal, or use `MvNormal(param, cv)` if `MvNormal(μ, D::Diagonal)` is supported).

---

### C. No check that read matrix size matches current model

**Where:** `fit.jl` → `get_propcv` reads a matrix from file; `make_structures` / `load_model` do not check that its size matches `length(fittedparam)`.

If you change the model (e.g. different `fittedparam` or number of parameters) between the run that wrote param-stats and the run that uses `propcv < 0`, the read covariance has the **wrong dimension**. Then `proposal_dist(param, cv, model)` either errors (dimension mismatch in `MvNormal`) or behaves in an undefined way.

**Suggestion:** When `propcv` is a tuple from file (or when it’s a matrix), in `make_structures` (after `set_fittedparam`) check `size(cv[1], 1) == length(fittedparam)` and throw a clear error if not. Optionally, in `get_propcv` you could return something that carries the dimension so the caller can validate even before having `fittedparam` (e.g. in a different entry point).

---

### D. Scaling may be too aggressive or inconsistent

**Where:** `fit.jl` → `get_propcv`: after reading the matrix, `cv = 2.38^2 * cv / size(cv, 1)`.

The factor `2.38^2/d` is the usual optimal scaling for a Gaussian random-walk proposal in dimension d. So the **formula** is correct. Possible issues:

- The **previous** run may have had different data, prior, or model, so the posterior covariance is not the same; reusing and scaling it is then heuristic.
- If the file covariance was already scaled (e.g. by some other code path), double-scaling could make proposals too small or too large.
- For hierarchical models, the same global scaling is applied even though only part of the parameters may use the matrix (hyperparameters); the rest use the scalar from the tuple. That might be fine, but it’s worth being aware of.

**Suggestion:** Document that the read matrix is treated as **unscaled** posterior covariance and is always rescaled by `2.38^2/n`. If any other code path ever writes a pre-scaled matrix, either stop doing that or add a flag so `get_propcv` does not scale again.

---

### E. Which block is read?

**Where:** `io.jl` → `read_bottom_float_block` returns the **last** contiguous block of numeric rows that form a square block (bottom `k` rows × first `k` columns).

`write_param_stats` writes, in order: labels, mean, std, med, mad, qparam, corparam, **covparam**, **covlogparam**. So the last block is **covlogparam** (covariance of **transformed** parameters). That matches the space of `param` used in the sampler. So this part is consistent; the only risk is if the file format or write order changes and the “bottom block” is no longer covlogparam.

---

### F. Annealing never updates the proposal

**Where:** `metropolis_hastings.jl` → `anneal` uses `proposalcv = model.proposal` and never updates `d` from a new covariance. So annealing always uses the same proposal (including your matrix). That’s consistent; the more impactful issue is warmup overwriting the matrix (see A).

---

## 3. Summary of suggested code changes

| Priority | Issue | Suggested change |
|----------|--------|-------------------|
| High | Warmup overwrites user covariance | Only adapt in warmup when initial proposal is scalar/vector; if it’s Matrix or Tuple, keep `model.proposal` for the main run (or add an option to “keep initial proposal”). |
| High | No `Diagonal` method | Add `proposal_dist(param::Vector, cv::Diagonal, model, indiv=0.001)` (e.g. delegate to `Matrix(cv)` or implement via product of Normals / `MvNormal`). |
| Medium | No dimension check | In `make_structures`, when proposal is Tuple or Matrix, check `size(cv, 1) == length(fittedparam)` and throw a clear error. |
| Low | Documentation | Document that file covariance is in transformed space, is always rescaled by `2.38^2/n`, and that changing the model between runs can make the read matrix invalid. |

---

## 4. Minimal code references

- **get_propcv (tuple from file):** `fit.jl` ~2947–2961  
- **Warmup overwriting:** `metropolis_hastings.jl` ~392–422  
- **proposal_dist(Matrix/Tuple):** `metropolis_hastings.jl` ~736–752  
- **Diagonal fallback:** `metropolis_hastings.jl` ~414–418  
- **Stats / covlogparam:** `metropolis_hastings.jl` ~939–962  
- **write_param_stats:** `io.jl` ~2772–2784  
- **read_bottom_float_block:** `io.jl` ~3955–3988  
