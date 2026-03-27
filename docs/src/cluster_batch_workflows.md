# Cluster and batch workflows

This chapter is for **everyone who installs StochasticGene.jl** (`Pkg.add("StochasticGene")`) and needs to:

- run many `fit` jobs on a cluster (**including NIH Biowulf**) using **swarm files** and [`makeswarm`](@ref);
- follow the **recommended coupled-model workflow**: fit **individual units first**, then **merge** those fitted rates into **one initial rate file** for the **coupled** model.

The implementations live in **`biowulf.jl`** (swarms, run-spec presets) and **`io.jl`** (merging rate tables). Function signatures and defaults are in the **docstrings**; this page is the **narrative** guide published with the [GitHub-hosted documentation](https://nih-niddk-mbs.github.io/StochasticGene.jl/stable/).

---

## Coupled models: recommended workflow (single units → merge → coupled fit)

For **coupled** transcriptional units (e.g. enhancer + gene, with or without a hidden unit), the **suggested workflow** is:

1. **Fit each unit as its own single-unit model**  
   Run separate fits (e.g. enhancer-only and gene-only traces or histograms) so each produces a standard MCMC **`rates_*.txt`** (rows = posterior samples, columns = rate headers for that unit).  
   Use normal [`fit`](@ref) calls or batch them with [`makeswarm_models`](@ref) / [`makeswarmfiles`](@ref) in **single-unit** mode.

2. **Merge the fitted rates into one wide table**  
   Stack the columns from the two (or more) unit files and append **coupling** placeholder columns using [`create_combined_file`](@ref) (two units) or [`create_combined_file_mult`](@ref) (more than two).  
   You choose **`Nenh` / `Ngene`** (or per-unit column counts) to match how each **set** of rates is laid out in your files (see docstrings).  
   For **many keys** (e.g. from a CSV of model names), use [`create_combined_files`](@ref) or [`create_combined_files_driver`](@ref), which call `create_combined_file` once per key and name outputs with [`combined_rates_key`](@ref).

3. **Run the coupled fit using the combined file as the starting rates**  
   Point **`infolder`** / **`inlabel`** (or your run spec) at that merged file so the coupled MCMC **warm-starts** from the stacked single-unit posteriors.  
   The coupled [`fit`](@ref) uses tuple **`G`**, **`R`**, **`coupling`**, joint datatype (e.g. **`tracejoint`**), etc. Coupling strengths are then **estimated** in the coupled run (the placeholder columns from step 2 get updated).

4. **Optional: batch everything on the cluster**  
   Use [`makeswarm`](@ref) or [`makeswarmfiles`](@ref) so each job runs `fit(; key=..., ...)` from prewritten **`info_<key>`** specs (see [Run specification (info TOML)](run_spec_toml.md)).

This order—**individual fits → merge → coupled fit**—is the standard way to get a sensible **initial** combined rate file for coupled models without fitting all parameters cold.

---

## NIH Biowulf: using [`makeswarm`](@ref)

[`makeswarm`](@ref) does **not** submit jobs to the scheduler by itself. It **writes files** you submit with Biowulf’s **`swarm`** (or your own `sbatch` wrappers):

- **`<swarmfile>.swarm`** — one command line per run key (each line runs `julia` with your project and one fit script).
- **`fitscript_<key>.jl`** per key — typically calls `fit(; key="<key>", ...)` with shared options (`resultfolder`, `maxtime`, `samplesteps`, etc.).

**Typical use on Biowulf**

1. Install StochasticGene in your Julia environment (see [Installation](installation.md), including [Biowulf Installation](installation.md#Biowulf-Installation-NIH-HPC)).
2. From Julia (interactive session or batch script), run something like:

```julia
using StochasticGene

makeswarm(
    ["runA", "runB"];           # keys; must match info_<key> / rates_<key> naming you use
    filedir      = "my_swarm",  # directory where .swarm and .jl files are written
    resultfolder = "my_results",
    root         = ".",
    project      = "/path/to/your/StochasticGene.jl",  # or "" if using the default environment
    nchains      = 4,
    nthreads     = 1,
    maxtime      = 72000.0,
    samplesteps  = 1_000_000,
)
```

3. Submit the swarm from the shell (example):

```bash
cd my_swarm
swarm -g 4 -t 16 -b 1 --time 24:00:00 --module julialang -f fit.swarm
```

Adjust **`-g`**, **`-t`**, **time**, and **module** to match your allocation and Julia module name on Biowulf.

**Generating keys and `info_<key>` in bulk**

- [`write_run_spec_preset`](@ref) — write `info_<key>.jld2` + marker TOML for one key.
- [`makeswarm_models`](@ref) — sweep single-unit `G,R,S,insertstep`, write presets, then call `makeswarm`.
- [`makeswarmfiles`](@ref) — **unified** entry: coupled key lists (CSV, explicit `base_keys`, or H3 grids) **or** single-unit sweeps; writes presets and runs `makeswarm`. See its docstring for the mutually exclusive modes.
- [`makeswarmfiles_h3_latent`](@ref) — convenience for H3 latent key grids.

---

## Key-based naming

Many batch helpers assume a string **`key`** per run:

- `results/<resultfolder>/info_<key>.toml` and `info_<key>.jld2`
- `rates_<key>.txt`

See [Run specification (info TOML)](run_spec_toml.md). Presets for cluster reruns are written with [`write_run_spec_preset`](@ref).

---

## Combined rate files (`io.jl`) — reference

- [`read_rates_table`](@ref), [`write_rates_table`](@ref)
- [`merge_coupled_two_unit_rates`](@ref), [`merge_coupled_stacked_units`](@ref)
- [`create_combined_file`](@ref), [`create_combined_file_mult`](@ref)
- [`read_combined_file_specs_csv`](@ref), [`create_combined_files_driver`](@ref), [`create_combined_files`](@ref), [`create_combined_files_h3_latent`](@ref)

---

## After the coupled fit

Post-processing examples: [`write_correlation_functions`](@ref), [`write_traces`](@ref), and other analysis functions in the [API Reference](api/index.md).

---

## See also

- [Run specification (info TOML)](run_spec_toml.md)
- [Coupled model analysis](examples/coupled_models.md) (example-focused; batch mechanics are on this page)
- [Installation](installation.md) (includes a **Biowulf** subsection)
- [Model fitting (`fit`)](api/fit.md)

---

## Maintainer note: where to document what

| Topic | Canonical place |
|-------|-------------------|
| User workflows, Biowulf, coupled merge order | **This page** (hosted docs) |
| `info_<key>` file format | [run_spec_toml.md](run_spec_toml.md) |
| README on GitHub | Short pointer + link to stable docs |
| Exact function signatures | Docstrings in `biowulf.jl` / `io.jl` |
