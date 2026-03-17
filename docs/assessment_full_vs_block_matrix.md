# Assessment: Full-matrix element expansion vs block-wise (RG) construction

## Element contract: summands, not final entries

Each `Element` does **not** represent the final matrix entry at (row,col). It is one **component** that must be summed over. We build each block/unit with its own element structure, then sum over all at assembly time. So: no products of rates, only sums; the element approach works by overloading "element" to mean "one term in the sum" for each (row,col).

## How unit-level elements work (overrepresentation)

At the unit level, each transition matrix entry is built as a **sum** of `Element` contributions:

- **`make_mat!(T, elements, rates)`**: `T[e.a, e.b] += e.pm * rates[e.index]` for each element.
- So each matrix entry is **sum over elements at (i,j) of (pm × rate[index])** — no products of rates.

You overrepresent by using **multiple elements that sum to one physical transition**:

1. **G block** (`set_elements_G!`): Each transition (e.g. 1→2) is two elements: diagonal `(1,1, γ, -1)` and off-diagonal `(2,1, γ, +1)`, so the net flow is a single rate γ.
2. **Replication into T** (`elements_TG!`, `elements_TR!`, `elements_RG!`): The same (index, pm) is replicated across blocks (e.g. G block repeated for each R block via `Element(e.a + j, e.b + j, e.index, e.pm)`).

So unit-level: **each (i,j) is a sum of rate terms only**. No products.

## What happens under the Kronecker sum (block-wise)

Uncoupled full matrix:

- **T_full = T1⊗I + I⊗T2** (and similarly for n units).
- So **T_full[(i,i'), (j,j')] = T1[i,j] δ(i',j') + δ(i,j) T2[i',j']**.
- Each entry is still a **sum of rates**: at most one term from unit 1 and one from unit 2. No products of rates from different units.

Coupling term:

- **coupling_strength[k] × (U_β ⊗ V_α)** — here the full entry is **coupling_strength × (U entry) × (V entry)**. V is built from target rates, so we do get a product (coupling_strength × rate). That is the only place with a product; we fixed it by storing coupling elements with target rate index and multiplying by coupling_strength at assembly.

## Full-matrix expansion in principle

Expansion only **replicates** unit elements into full (row, col):

- For each unit-α element (a, b, index, pm), we generate one full-space element per state of the other slots: (row, col, rate_base + index - 1, pm).
- So each full (row, col) still gets a **sum** of (pm × rate) terms — the same as in the block-wise case. We do **not** introduce products of rates by expanding.

So in principle the full-matrix representation is still “elements are sums of rates”; the bug is almost certainly in **how** we set (row, col) or rate index during expansion (e.g. indexing, ordering, or double-counting), not in the fact that we use a flat element list.

## Where complexity does explode

1. **Number of elements**: We go from O(nT_α) elements per unit to O(N × n_elems_α) with N = ∏ nT. So the list gets big and indexing is easy to get wrong.
2. **Coupling**: Coupling entries are the only ones that are products (strength × source × target rate). We handled that by a separate `elements_coupling` list and a second pass in assembly.
3. **Debugging**: With a single flat list, tracing “which unit/slot contributed this (row, col)?” is harder than in the block-wise view where T1, T2, and Kronecker steps are explicit.

## Options

### A. Keep block-wise (RG) method

- **Pros**: Already correct for 2 units; structure is clear (T1, T2, kron_left/kron_right); no expansion bugs.
- **Cons**: 3-unit path currently hits DimensionMismatch in legacy `make_mat_TC`; coupling is limited to “one block at a time” (hard to do R-any, combined sources/targets, etc.).

**If we keep block-wise**: Fix the 3-unit DimensionMismatch in the RG path and use it for all n-unit uncoupled + current coupling. No full-matrix expansion for the uncoupled part.

### B. Fix the full-matrix expansion

- **Pros**: Single list of elements; one assembly loop; natural for “scan full space” coupling (R-any, combined, etc.).
- **Cons**: Need to find and fix the indexing bug (e.g. (73,73) wrong by 1.325); expansion logic is subtle (slot order, other_k → state mapping).

**If we fix expansion**: Add targeted tests (e.g. compare one row/column of T_full vs T_RG, or count contributions per (row,col)) to pin down the wrong (row, col) or rate index.

### C. Hybrid

- **Uncoupled**: Always use block-wise (build T[α], then Tc = ∑_α kron_expand(T[α])). No flat uncoupled element list.
- **Coupling**: Use full-space coupling elements only (scan full state space and push (row, col, rate_index, pm) for each coupling connection). Assembly: uncoupled from block-wise, then add coupling from elements.

**Pros**: Uncoupled part is proven and simple; coupling can still be general (R-any, etc.). **Cons**: Two code paths (block-wise for T, element list for coupling); need to combine them in one T matrix.

## Recommendation

- **Short term**: Treat the full-matrix **uncoupled** expansion as the place where the bug is (element setting), not the representation itself. So: either **fix the expansion** (Option B) with focused tests, or **drop uncoupled expansion** and keep **block-wise for uncoupled** (Option A or C).
- **If we want R-any / complex coupling soon**: Prefer **Option C** (block-wise uncoupled + full-space coupling only), then we don’t depend on fixing the uncoupled expansion.
- **If we want a single, simple “one list, one loop” design**: Then **Option B** (fix expansion) is the right target; the assessment above says the design is sound and the issue is in the indexing/replication in `expand_unit_elements_to_full` (or in how we build the legacy comp used there).
