# Cluster and batch workflows

This chapter is for **everyone who installs StochasticGene.jl** (`Pkg.add("StochasticGene")`) and needs to:

- run many `fit` jobs on a cluster (**including NIH Biowulf**) using **swarm files** and the composable helpers [`emit_fitscript`](@ref) / [`emit_fitscripts`](@ref) + [`emit_swarm_batch`](@ref) (legacy [`makeswarm`](@ref) / [`makeswarmfiles`](@ref) are deprecated but still work);
- follow the **recommended coupled-model workflow**: fit **individual units first**, then **merge** those fitted rates into **one initial rate file** for the **coupled** model.

The implementations live in **`biowulf.jl`** (swarms, run-spec presets) and **`io.jl`** (merging rate tables). Function signatures and defaults are in the **docstrings**; this page is the **narrative** guide published with the [GitHub-hosted documentation](https://nih-niddk-mbs.github.io/StochasticGene.jl/stable/).

## Composable API (preferred)

Use these in new code; legacy names remain for compatibility but emit `Base.depwarn`:

| Step | Functions |
|------|-----------|
| Write fit scripts | [`emit_fitscript`](@ref), [`emit_fitscripts`](@ref) (set `write_swarm=true` to also emit a `.swarm`) |
| Write swarm files only | [`emit_swarm`](@ref), [`emit_swarm_batch`](@ref) (`per_key_scripts=true` matches “one `.jl` per key”) |
| Stage / migrate folders | [`stage_run`](@ref) (alias [`setup`](@ref)): `specs/`, `scripts/`, `swarm/`, `inputs/`, optional [`create_combined_file`](@ref) via `combined=` |
| Defaults / profiles | [`load_stochasticgene_toml_layers`](@ref) reads `~/.stochasticgene.toml` and `root/.stochasticgene.toml`, then optional `[profiles.NAME]` |
| Model grids / specs | [`model_grid`](@ref), [`fit_specs`](@ref), [`default_priors`](@ref), [`coupling_grid`](@ref) |

**Gene panels:** [`write_fitfile_genes`](@ref) and [`makeswarm_genes`](@ref) produce **one** shared script with `gene=ARGS[1]` and one swarm line per gene (still the recommended pattern for many genes, one model).

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
   Point **`datapath`** / **`resultfolder`** / **`label`** (or your run spec) at that merged file so the coupled MCMC **warm-starts** from the stacked single-unit posteriors.  
   The coupled [`fit`](@ref) uses tuple **`G`**, **`R`**, **`coupling`**, joint datatype (e.g. **`tracejoint`**), etc. Coupling strengths are then **estimated** in the coupled run (the placeholder columns from step 2 get updated).

4. **Optional: batch everything on the cluster**  
   Use [`emit_fitscripts`](@ref) (with `write_swarm=true`) or legacy [`makeswarm`](@ref) / [`makeswarmfiles`](@ref) so each job runs `fit(; key=..., ...)` from prewritten **`info_<key>`** specs (see [Run specification (info TOML)](run_spec_toml.md)).

This order—**individual fits → merge → coupled fit**—is the standard way to get a sensible **initial** combined rate file for coupled models without fitting all parameters cold.

---

## NIH Biowulf: key-based jobs ([`emit_fitscripts`](@ref) / legacy [`makeswarm`](@ref))

[`emit_fitscripts`](@ref) and legacy [`makeswarm`](@ref) **do not** submit jobs to the scheduler by themselves. They **write files** you submit with Biowulf’s **`swarm`** (or your own `sbatch` wrappers):

- **`<swarmfile>.swarm`** — one command line per run key (each line runs `julia` with your project and one fit script).
- **`fitscript_<key>.jl`** per key — typically calls `fit(; key="<key>", ...)` with shared options (`resultfolder`, `maxtime`, `samplesteps`, `warmupsteps`, `inference_method`, `device`, `parallel`, `gradient`, etc.). Symbol-valued overrides (e.g. `inference_method=:nuts`) are supported in minimal scripts; full key-based specs store the same keys in **`info_<key>.jld2`** (see [Run specification](run_spec_toml.md)).

**Typical use on Biowulf**

1. Install StochasticGene in your Julia environment (see [Installation](installation.md), including [Biowulf Installation](installation.md#Biowulf-Installation-NIH-HPC)).
2. From Julia (interactive session or batch script), run something like:

```julia
using StochasticGene

# Legacy one-liner (still works; emits depwarn):
makeswarm(
    ["runA", "runB"];
    filedir      = "my_swarm",
    resultfolder = "my_results",
    root         = ".",
    project      = "/path/to/your/StochasticGene.jl",
    nchains      = 4,
    nthreads     = 1,
    maxtime      = 72000.0,
    samplesteps  = 1_000_000,
)

# Preferred split: write scripts with `emit_fitscript(path, key; ...)` then
# `emit_swarm_batch(...; per_key_scripts=true)` — see docstrings.
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

### Swarm `julia -p`, `nchains`, merged `info_<key>`, and `root`

- **Parallel workers:** The swarm command should use **`-p N`** (or equivalent) consistent with how many chains run in parallel. For [`makeswarmfiles`](@ref) / [`makeswarm_models`](@ref), if you do **not** pass an explicit swarm-only **`nchains=`** in kwargs, the generated **`-p`** is taken from each run spec’s **`nchains`** (e.g. coupled defaults often use 16), so it stays aligned with **`fit(; …, nchains=…)`**. See the [`makeswarmfiles`](@ref) docstring. For **NUTS/ADVI**, `nchains` still controls how many independent chains are launched; within-chain parallelism follows each method’s options (`parallel` on [`NUTSOptions`](@ref) / [`ADVIOptions`](@ref), set via [`load_options`](@ref) from the run spec).

- **Merged presets:** With **`merge_existing_info=true`** (default), older **`info_<key>.jld2`** files are merged into new specs. Legacy **`trace_specs`** sometimes used a huge **`t_end`** (historical “open end” sentinel). When saving, [`write_run_spec_preset`](@ref) runs [`normalize_trace_specs_legacy_t_end!`](@ref) so those values are rewritten to **`t_end = -1.0`**, matching current [`default_trace_specs_for_coupled`](@ref) and avoiding invalid frame indices in [`read_tracefiles`](@ref).

- **`root` in generated fit scripts:** Scripts list **`root`** exactly as in the run spec (no forced `abspath`). Use **`root="."`** if the job’s **working directory** is the project root (set **`cd`** in the swarm or submit from the right folder). Paths resolved in an **interactive** Biowulf session can differ from **batch** jobs; **`"."`** avoids baking in an interactive-only absolute path.

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
