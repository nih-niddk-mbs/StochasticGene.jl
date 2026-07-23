# StochasticGene.jl

**Version 1.11.0**

[![Documents](https://img.shields.io/badge/Documents-blue.svg)](https://nih-niddk-mbs.github.io/StochasticGene.jl/dev/)

Julia package for **stochastic models of gene transcription**, **Bayesian inference** (Metropolis–Hastings, NUTS, and mean-field ADVI on the same transformed parameters), and **analysis** of smFISH, scRNA-seq, live-cell traces, and dwell-time data.

---

## Documentation

| Resource | URL |
|----------|-----|
| **Manual (hosted)** | [Documents](https://nih-niddk-mbs.github.io/StochasticGene.jl/dev/) |
| **Cluster & batch workflows** | [Cluster and batch workflows](https://nih-niddk-mbs.github.io/StochasticGene.jl/dev/cluster_batch_workflows.html) — Biowulf swarms, `makeswarm`, `makeswarmfiles`, combined rate files, key-based `fit` |
| **Package overview** | [Repository layout, nomenclature, data types](https://nih-niddk-mbs.github.io/StochasticGene.jl/dev/package_overview.html) |
| **Issues** | [GitHub Issues](https://github.com/nih-niddk-mbs/StochasticGene.jl/issues) |

---

## Installation

From the Julia REPL, enter package mode with `]` and run:

```julia
pkg> add StochasticGene
```

Julia **≥ 1.11** is required for the v1.11 beta / 2.0-beta line (see `Project.toml`). For development, clone the repo and `pkg> dev path/to/StochasticGene.jl`.

---

## What this package provides

- **Models:** Generalized telegraph / GRSM dynamics — arbitrary **G** (gene) states, **R** pre-RNA steps, **S** splice sites, reporter **insertstep**, multiple alleles, **coupled** transcribing units (e.g. enhancer–gene, hidden latent units).
- **Inference:** `fit` selects **MH**, **NUTS**, or **ADVI** via `inference_method`; shared budgets (`samplesteps`, `warmupsteps`, …) map through `load_options` to method-specific structs. MH uses fixed proposal CVs or saved covariance proposals and multi-chain `Distributed` pooling; cluster helpers emit compatible `fit(; …)` overrides (see **Cluster & batch workflows** in the manual).
- **Data:** Stationary RNA histograms, ON/OFF and dwell-time histograms, intensity traces (`.trk` etc.), joint traces, grids — alone or in combination. In the v1.11 beta / 2.0-beta line, multimodal fits use `CombinedData` via tuple/vector `datatype` values such as `(:rna, :dwelltime)`.
- **Batch & HPC:** Helpers in `biowulf.jl` write **swarm** files and **fit scripts** for NIH Biowulf or any scheduler; run specs can be saved as `info_<key>.toml` + `info_<key>.jld2` for reproducible `fit(; key=...)`.

---

## Repository layout (typical project)

When you work in a **project root** (the `root` argument to `fit`, often `"."`):

- **`data/`** — inputs: per-gene files, trace folders, `*_alleles.csv`, `*_halflife.csv`, etc.
- **`results/<resultfolder>/`** — outputs: rates, measures, `info_<key>.*`, etc.

The **`results/`** directory is listed in **`.gitignore`** so local fits are not committed. Keep raw outputs out of git; share summaries, figures, or exported tables instead.

Helper: `rna_setup("myproject")` creates a sensible `data/` / `results/` skeleton (see the manual).

---

## Quick start (RNA histogram)

For `datatype="rna"`, data are steady-state RNA count histograms stored as one
text file per gene/condition:

```text
data/HCT116_testdata/
├── CENPL_MOCK.txt
└── MYC_MOCK.txt
```

Each file's first column is the histogram over RNA copy number bins.

```julia
using StochasticGene

fits, stats, measures, data, model, options = fit(
    G = 2,
    R = 0,
    S = 0,
    insertstep = 1,
    transitions = ([1, 2], [2, 1]),
    datatype = "rna",
    datapath = "data/HCT116_testdata",
    gene = "MYC",
    cell = "HCT116",
    datacond = "MOCK",
    resultfolder = "HCT116_test",
)
```

Defaults in `fit()` point at bundled mock data if paths match the package examples.

For multimodal data, prefer the `CombinedData` API:

```julia
fits, stats, measures, data, model, options = fit(
    datatype = (:rna, :dwelltime),
    datapath = (
        rna = "HCT116_smFISH",
        dwelltime = ["dwell/MYC_ON.csv", "dwell/MYC_OFF.csv"],
    ),
    gene = "MYC",
    cell = "HCT116",
    datacond = "MOCK",
    resultfolder = "HCT116_test",
    dwell_specs = [(unit = 1, onstates = [[2], Int[]], dttype = ["ON", "OFF"])],
)
```

Legacy `infolder` / `inlabel` style inputs have been retired. Use `root`, `datapath`, `label`, and `resultfolder` directly; `datapath` entries are resolved under `root/data` when needed, and outputs are written under `root/results/<resultfolder>`.

---

## Key-based runs and batch jobs

For reproducible configs, save a run spec and call **`fit(; key="33il", resultfolder="my_results", root=".")`**. The package looks for **`results/.../info_<key>.toml`** (with JLD2 companion) and merges with keyword overrides.

Batch helpers:

- **`makeswarm(["key1", "key2"]; filedir=..., resultfolder=..., ...)`** — one swarm line and **`fitscript_<key>.jl`** per key (`fit(; key=..., ...)`).
- **`makeswarm_genes(["GENE1", "GENE2"]; ...)`** — same model, one job per **gene** (genome-scale scRNA-style).
- **`makeswarm_genes(; datapath="data/HCT116_testdata", datacond="MOCK", ...)`** — scan the RNA data folder for matching `GENE_COND.txt` files, then run the same gene-panel writer. This is the current replacement for the old v0.7-style `makeswarm(; datafolder=..., conds=...)` RNA sweep. By default it uses `checkgenes` and restricts to genes with available halflife/allele metadata; pass `filter_metadata=false` to scan filenames only.
- **`makeswarm_models`** / **`makeswarmfiles`** — model sweeps, coupled CSV workflows, combined-rate keys; see the [cluster & batch chapter](https://nih-niddk-mbs.github.io/StochasticGene.jl/stable/cluster_batch_workflows.html).

The `makeswarm` family always writes a plain command list (`fit.swarm` by default).
Use `scheduler` to add a launcher for the machine you are on:

```julia
# Biowulf default: submit fit.swarm with swarm
makeswarm_genes(; datapath="data/HCT116_testdata", datacond="MOCK")

# Slurm: writes fit.swarm and fit_slurm.sh
makeswarm_genes(; datapath="data/HCT116_testdata", datacond="MOCK",
    scheduler=:slurm, scheduler_jobs=100, slurm_mem="8G", slurm_time="02:00:00")

# GNU Parallel: writes fit.swarm and fit_parallel.sh
makeswarm_genes(; datapath="data/HCT116_testdata", datacond="MOCK",
    scheduler=:parallel, scheduler_jobs=16)
```

Then run `swarm -f fit.swarm` on Biowulf, `sbatch fit_slurm.sh` on Slurm,
`./fit_parallel.sh` with GNU Parallel, or `bash fit.swarm` for a sequential local run.

**Coupled models:** fit each unit separately, merge rates with **`create_combined_file`** / **`create_combined_file_mult`**, then run the coupled fit warm-started from that table (detailed steps in the manual).

---

## Source code map (`src/`)

| File | Role |
|------|------|
| `common.jl` | Types and shared utilities |
| `transition_rate_*.jl` | Rate matrices and master-equation structure |
| `likelihoods.jl`, `chemical_master.jl` | Likelihoods and CME backends |
| `metropolis_hastings.jl` | MCMC |
| `fit.jl` | **`fit`**, run specs, defaults |
| `io.jl` | I/O, combined files, run-spec read/write |
| `utilities.jl` | Helpers |
| `simulator_coupled.jl` | **`simulator`**, stochastic simulation |
| `analysis.jl` | Post-fit analysis, correlation functions, exports |
| `hmm.jl` | Observation models / HMM pieces |
| `coupled_csv.jl` | **`Coupled_models_to_test.csv`** → coupling specs for batch coupled fits |
| `biowulf.jl` | **`makeswarm`**, **`makeswarmfiles`**, swarm writers |
| `test.jl` | Regression / dev tests |

Public API is exported from `StochasticGene.jl`; full lists appear under **API** in the hosted manual.

---

## Units and nomenclature

Rates are in **inverse minutes**; halflife files use **hours** where applicable. The package distinguishes **gene states** (G) from **R steps** along the gene (occupancy patterns are combinatorial — see [Package overview](https://nih-niddk-mbs.github.io/StochasticGene.jl/stable/package_overview.html)).

---

## NIH Biowulf

Load Julia (`module load julialang`), generate swarm files from the REPL with **`makeswarm`** / **`makeswarm_genes`**, then submit with **`swarm`** (allocate time ≥ `maxtime`, enough `-t` / memory for `nchains`). If **`method`** is passed through swarm-generating functions, it must be a **String** in the generated script (e.g. `"Tsit5()"`) so the file parses correctly — see docstrings in `biowulf.jl`.

Example RNA folder sweep:

```julia
using StochasticGene

makeswarm_genes(
    datapath="data/HCT116_testdata",
    datacond="MOCK",
    resultfolder="HCT116_test",
    filedir="run-HCT116-testdata-rna",
    cell="HCT116",
    datatype="rna",
    G=2,
    R=0,
    S=0,
    insertstep=1,
    transitions=([1, 2], [2, 1]),
    nchains=8,
    project="/home/carsonc/github/StochasticGene.jl/",
)
```

For genome-scale scRNA folders, `makeswarm_genes` now writes one command file by
default, with one line per gene. On Biowulf, keep that single file and use
`swarm -b` at submission time to bundle many gene commands into fewer Slurm
array tasks:

```julia
out = makeswarm_genes(
    root="/data/carsonc/scrna",
    datapath="RamosNELFA_NEG_IFNa_rep1_Sdata",
    datacond="NIr1",
    cell="U3A",
    resultfolder="RamosNELFA_NEG_IFNa_rep1",
    filedir="run-RamosNELFA_NEG_IFNa_rep1",
    nchains=2,
    nthreads=1,
    project="/home/carsonc/github/StochasticGene.jl/",
)
```

Submit on Biowulf with bundling, choosing `-b` so the bundled walltime remains
reasonable:

```bash
cd run-RamosNELFA_NEG_IFNa_rep1
swarm -f fit.swarm -b 20 -g 4 -t 2 --time 2:00:00 --merge-output --module julia
```

Set `batchsize=<N>` only if you deliberately want multiple command files.

On non-Biowulf systems, pass `scheduler=:slurm` to also write `fit_slurm.sh`,
or `scheduler=:parallel` to also write `fit_parallel.sh`.

After the jobs finish, assemble result tables:

```julia
write_dataframes_only(
    "results/HCT116_test",
    "data/HCT116_testdata";
    datatype="rna",
)
```

For key-based folders containing `rates_<key>.txt`, `param-stats_<key>.txt`,
and `info_<key>.jld2`, write one key summary with:

```julia
write_dataframes_only_key("results/HCT116_test"; datatype="rna")
```

---

## Citation

If you use this package, please cite the relevant papers for your application, for example:

- Rodriguez *et al.*, *Cell* (2018)  
- Wan *et al.*, *Cell* (2021)  
- Trzaskoma *et al.*, *Science Advances* (2024)

---

## License

MIT License — see [LICENSE](LICENSE).

---

## Contributing

See [Contributing](https://nih-niddk-mbs.github.io/StochasticGene.jl/dev/contributing.html) in the manual. Pull requests and issue reports are welcome.
