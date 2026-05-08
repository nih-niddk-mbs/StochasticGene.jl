# Cluster and batch workflows

This chapter is for **everyone who installs StochasticGene.jl** (`Pkg.add("StochasticGene")`) and needs to:

- run many `fit` jobs on a cluster (**including NIH Biowulf**) using **command files** (scheduler-agnostic) or **swarm files** (Biowulf naming) from the stage helpers [`make_fitscript`][make_fitscript], [`make_fitscripts_from_csv`][make_fitscripts_from_csv], [`make_commandfile`][make_commandfile], [`make_commandfile_from_csv`][make_commandfile_from_csv], [`make_fitscripts_and_commandfile_from_csv`][make_fitscripts_and_commandfile_from_csv], and compatibility wrappers [`make_swarmfile_from_csv`][make_swarmfile_from_csv], [`make_fitscripts_and_swarm_from_csv`][make_fitscripts_and_swarm_from_csv];
- follow the **recommended coupled-model workflow**: fit **individual units first**, then **merge** those fitted rates into **one initial rate file** for the **coupled** model.

The implementations live in **`stage.jl`** (stage-native scripts + command files), **`biowulf.jl`** (legacy and Biowulf-specific wrappers), and **`io.jl`** (merging rate tables). Function signatures and defaults are in the **docstrings**; this page is the **narrative** guide published with the [GitHub-hosted documentation](https://nih-niddk-mbs.github.io/StochasticGene.jl/dev/).

Source entry points:
- [`src/stage.jl`](https://github.com/nih-niddk-mbs/StochasticGene.jl/blob/main/src/stage.jl)
- [`src/biowulf.jl`](https://github.com/nih-niddk-mbs/StochasticGene.jl/blob/main/src/biowulf.jl)
- [`src/io.jl`](https://github.com/nih-niddk-mbs/StochasticGene.jl/blob/main/src/io.jl)

## Stage API (preferred)

Use these in new code:

| Step | Functions |
|------|-----------|
| Write one key script | [`make_fitscript`][make_fitscript] |
| Write scripts from keys CSV | [`make_fitscripts_from_csv`][make_fitscripts_from_csv] |
| Build one launch command string | [`build_julia_script_command`][build_julia_script_command] |
| Write command file from script list | [`make_commandfile`][make_commandfile], [`write_julia_command_file`][write_julia_command_file] |
| Write command file from keys CSV | [`make_commandfile_from_csv`][make_commandfile_from_csv] |
| Write scripts + command file from CSV | [`make_fitscripts_and_commandfile_from_csv`][make_fitscripts_and_commandfile_from_csv] |
| Biowulf naming wrappers | [`make_swarmfile_from_csv`][make_swarmfile_from_csv], [`make_fitscripts_and_swarm_from_csv`][make_fitscripts_and_swarm_from_csv] |

**Gene panels:** [`write_fitfile_genes`][write_fitfile_genes] and [`makeswarm_genes`][makeswarm_genes] produce **one** shared script with `gene=ARGS[1]` and one swarm line per gene (still the recommended pattern for many genes, one model).

---

## Coupled models: recommended workflow (single units → merge → coupled fit)

For **coupled** transcriptional units (e.g. enhancer + gene, with or without a hidden unit), the **suggested workflow** is:

1. **Fit each unit as its own single-unit model**  
   Run separate fits (e.g. enhancer-only and gene-only traces or histograms) so each produces a standard MCMC **`rates_*.txt`** (rows = posterior samples, columns = rate headers for that unit).  
   Use normal [`fit`][fit] calls or batch them with [`makeswarm_models`][makeswarm_models] / [`makeswarmfiles`][makeswarmfiles] in **single-unit** mode.

2. **Merge the fitted rates into one wide table**  
   Stack the columns from the two (or more) unit files and append **coupling** placeholder columns using [`create_combined_file`][create_combined_file] for two-unit models, or [`create_combined_file_mult`][create_combined_file_mult] for models with more than two units.  
   You choose **`Nenh` / `Ngene`** (or per-unit column counts) to match how each **set** of rates is laid out in your files (see docstrings).  
   For **many keys** (e.g. from a CSV of model names), use [`create_combined_files`][create_combined_files] or [`create_combined_files_driver`][create_combined_files_driver], which call `create_combined_file` once per key and name outputs with [`combined_rates_key`][combined_rates_key].

3. **Run the coupled fit using the combined file as the starting rates**  
   Point **`datapath`** / **`resultfolder`** / **`label`** (or your run spec) at that merged file so the coupled MCMC **warm-starts** from the stacked single-unit posteriors.  
   The coupled [`fit`][fit] uses tuple **`G`**, **`R`**, **`coupling`**, joint datatype (e.g. **`tracejoint`**), etc. Coupling strengths are then **estimated** in the coupled run (the placeholder columns from step 2 get updated).

4. **Optional: batch everything on the cluster**  
   Use stage-native [`make_fitscripts_from_csv`][make_fitscripts_from_csv] + [`make_commandfile_from_csv`][make_commandfile_from_csv] (or [`make_fitscripts_and_commandfile_from_csv`][make_fitscripts_and_commandfile_from_csv]) so each job runs `fit(; key=..., ...)` from prewritten **`info_<key>`** specs (see [Run specification (info TOML)](run_spec_toml.md)). For Biowulf-style `*.swarm` naming, use [`make_swarmfile_from_csv`][make_swarmfile_from_csv].

This order—**individual fits → merge → coupled fit**—is the standard way to get a sensible **initial** combined rate file for coupled models without fitting all parameters cold.

---

## Scheduler launch files and Biowulf swarm

Stage helpers **do not** submit jobs to a scheduler by themselves. They **write files** you submit with Biowulf’s **`swarm`** or another launcher (`sbatch`, GNU parallel, custom wrappers):

- **`<commandfile>`** (default `fit.commands`) — one command line per run key (`julia ... fitscript_<key>.jl`).
- **`<swarmfile>.swarm`** — same content style, Biowulf-compatible extension via [`make_swarmfile_from_csv`][make_swarmfile_from_csv].
- **`fitscript_<key>.jl`** per key — typically calls `fit(; key="<key>", ...)` with shared options (`resultfolder`, `maxtime`, `samplesteps`, `warmupsteps`, `inference_method`, `device`, `parallel`, `gradient`, etc.).

**Typical use on Biowulf**

1. Install StochasticGene in your Julia environment (see [Installation](installation.md), including [Biowulf Installation](installation.md#Biowulf-Installation-NIH-HPC)).
2. From Julia (interactive session or batch script), run something like:

```julia
using StochasticGene

# Stage-native pipeline: scripts + command file from a keys CSV.
# CSV requires a key column, default: Model_name
make_fitscripts_and_commandfile_from_csv(
    "keys.csv";
    filedir      = "my_jobs",
    commandfile  = "fit.commands",
    juliafile    = "fitscript",
    project      = "/path/to/your/StochasticGene.jl",
    nprocs       = 4,
    nthreads     = 1,
    resultfolder = "my_results",
    root         = ".",
    maxtime      = 72000.0,
    samplesteps  = 1_000_000,
)

# If you specifically want .swarm naming:
# make_fitscripts_and_swarm_from_csv("keys.csv"; filedir="my_jobs", swarmfile="fit", ...)
```

3. Submit the swarm from the shell (example):

```bash
cd my_jobs
swarm -g 4 -t 16 -b 1 --time 24:00:00 --module julialang -f fit.swarm
```

Adjust **`-g`**, **`-t`**, **time**, and **module** to match your allocation and Julia module name on Biowulf.

**Generating keys and `info_<key>` in bulk**

- [`write_run_spec_preset`][write_run_spec_preset] — write `info_<key>.jld2` + marker TOML for one key.
- [`makeswarm_models`][makeswarm_models] — sweep single-unit `G,R,S,insertstep`, write presets, then call Biowulf-oriented writers.
- [`makeswarmfiles`][makeswarmfiles] — **unified** legacy Biowulf entry: coupled key lists (CSV, explicit `base_keys`, or H3 grids) **or** single-unit sweeps; writes presets and emits swarm+scripts.
- [`makeswarmfiles_h3_latent`][makeswarmfiles_h3_latent] — convenience for H3 latent key grids.

### Swarm `julia -p`, `nchains`, merged `info_<key>`, and `root`

- **Parallel workers:** The swarm command should use **`-p N`** (or equivalent) consistent with how many chains run in parallel. For [`makeswarmfiles`][makeswarmfiles] / [`makeswarm_models`][makeswarm_models], if you do **not** pass an explicit swarm-only **`nchains=`** in kwargs, the generated **`-p`** is taken from each run spec’s **`nchains`** (e.g. coupled defaults often use 16), so it stays aligned with **`fit(; …, nchains=…)`**. See the [`makeswarmfiles`][makeswarmfiles] docstring. For **NUTS/ADVI**, `nchains` still controls how many independent chains are launched; within-chain parallelism follows each method’s options (`parallel` on [`NUTSOptions`][NUTSOptions] / [`ADVIOptions`][ADVIOptions], set via [`load_options`][load_options] from the run spec).

- **Merged presets:** With **`merge_existing_info=true`** (default), older **`info_<key>.jld2`** files are merged into new specs. Legacy **`trace_specs`** sometimes used a huge **`t_end`** (historical “open end” sentinel). When saving, [`write_run_spec_preset`][write_run_spec_preset] runs [`normalize_trace_specs_legacy_t_end!`][normalize_trace_specs_legacy_t_end] so those values are rewritten to **`t_end = -1.0`**, matching current [`default_trace_specs_for_coupled`][default_trace_specs_for_coupled] and avoiding invalid frame indices in [`read_tracefiles`][read_tracefiles].

- **`root` in generated fit scripts:** Scripts list **`root`** exactly as in the run spec (no forced `abspath`). Use **`root="."`** if the job’s **working directory** is the project root (set **`cd`** in the swarm or submit from the right folder). Paths resolved in an **interactive** Biowulf session can differ from **batch** jobs; **`"."`** avoids baking in an interactive-only absolute path.

---

## Key-based naming

Many batch helpers assume a string **`key`** per run:

- `results/<resultfolder>/info_<key>.toml` and `info_<key>.jld2`
- `rates_<key>.txt`

See [Run specification (info TOML)](run_spec_toml.md). Presets for cluster reruns are written with [`write_run_spec_preset`][write_run_spec_preset].

---

## Combined rate files (`io.jl`) — reference

- [`read_rates_table`][read_rates_table], [`write_rates_table`][write_rates_table]
- [`merge_coupled_two_unit_rates`][merge_coupled_two_unit_rates], [`merge_coupled_stacked_units`][merge_coupled_stacked_units]
- [`create_combined_file`][create_combined_file], [`create_combined_file_mult`][create_combined_file_mult]
- [`read_combined_file_specs_csv`][read_combined_file_specs_csv], [`create_combined_files_driver`][create_combined_files_driver], [`create_combined_files`][create_combined_files], [`create_combined_files_h3_latent`][create_combined_files_h3_latent]

---

## After the coupled fit

Post-processing examples: [`write_correlation_functions`][write_correlation_functions], [`write_traces`][write_traces], and other analysis functions in the [API Reference](api/index.md).

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
| User workflows, stage command files, Biowulf submission, coupled merge order | **This page** (hosted docs) |
| `info_<key>` file format | [run_spec_toml.md](run_spec_toml.md) |
| README on GitHub | Short pointer + link to hosted documents |
| Exact function signatures | Docstrings in `stage.jl` / `biowulf.jl` / `io.jl` |

[ADVIOptions]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=ADVIOptions&type=code
[NUTSOptions]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=NUTSOptions&type=code
[build_julia_script_command]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=build_julia_script_command&type=code
[combined_rates_key]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=combined_rates_key&type=code
[create_combined_file]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=create_combined_file&type=code
[create_combined_file_mult]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=create_combined_file_mult&type=code
[create_combined_files]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=create_combined_files&type=code
[create_combined_files_driver]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=create_combined_files_driver&type=code
[create_combined_files_h3_latent]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=create_combined_files_h3_latent&type=code
[default_trace_specs_for_coupled]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=default_trace_specs_for_coupled&type=code
[fit]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=fit&type=code
[load_options]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=load_options&type=code
[make_commandfile]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=make_commandfile&type=code
[make_commandfile_from_csv]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=make_commandfile_from_csv&type=code
[make_fitscript]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=make_fitscript&type=code
[make_fitscripts_and_commandfile_from_csv]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=make_fitscripts_and_commandfile_from_csv&type=code
[make_fitscripts_and_swarm_from_csv]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=make_fitscripts_and_swarm_from_csv&type=code
[make_fitscripts_from_csv]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=make_fitscripts_from_csv&type=code
[make_swarmfile_from_csv]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=make_swarmfile_from_csv&type=code
[makeswarm_genes]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=makeswarm_genes&type=code
[makeswarm_models]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=makeswarm_models&type=code
[makeswarmfiles]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=makeswarmfiles&type=code
[makeswarmfiles_h3_latent]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=makeswarmfiles_h3_latent&type=code
[merge_coupled_stacked_units]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=merge_coupled_stacked_units&type=code
[merge_coupled_two_unit_rates]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=merge_coupled_two_unit_rates&type=code
[normalize_trace_specs_legacy_t_end]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=normalize_trace_specs_legacy_t_end%21&type=code
[read_combined_file_specs_csv]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=read_combined_file_specs_csv&type=code
[read_rates_table]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=read_rates_table&type=code
[read_tracefiles]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=read_tracefiles&type=code
[write_correlation_functions]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=write_correlation_functions&type=code
[write_fitfile_genes]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=write_fitfile_genes&type=code
[write_julia_command_file]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=write_julia_command_file&type=code
[write_rates_table]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=write_rates_table&type=code
[write_run_spec_preset]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=write_run_spec_preset&type=code
[write_traces]: https://github.com/nih-niddk-mbs/StochasticGene.jl/search?q=write_traces&type=code
