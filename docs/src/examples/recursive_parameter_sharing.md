# Recursive Parameter-Sharing Trace Fits

Recursive parameter sharing lets one fit combine trace datasets that live in different folders while sharing some parameter families and splitting others. This is useful when the same biological model is fit to several related samples, for example 3Prime and 5Prime traces, CRISPR perturbations, or replicate groups.

This is not a hierarchical Bayesian shrinkage model: no learned hyperparameters are introduced. Group-specific or trace-specific copies are independent fitted parameters under the ordinary prior template. The fit keyword is still named `hierarchical` for compatibility with the existing `fit` API, but the object passed there is a recursive parameter-sharing specification.

The recursive interface is keyed: each dataset row carries metadata, the grouping tree reads those metadata keys, and `parameter_scope` decides which level owns each family of parameters.

## Minimal Example

For a simple two-sample fit with shared `G` transitions and initiation, but group-specific `R` transitions and group-specific noise:

```text
top
  group
```

```julia
using StochasticGene

grouping = HierarchySpec(
    HierarchyNode(:top, children=[
        HierarchyNode(:group, key=:sample, levels=[:ThreePrime, :FivePrime])
    ]);
    parameter_scope=Dict(
        :G => :top,
        :initiation => :top,
        :R => :group,
        :decay => :top,
        :noise => :group,
    ),
)
```

Each dataset declares its group through `metadata.sample`:

```julia
datasets = DatasetSpec[
    DatasetSpec(
        :ThreePrime_gene,
        "trace",
        "data/3Prime_gene_enhancer/including_background/short",
        "gene";
        metadata=(sample=:ThreePrime,),
        trace_specs=[(unit=1, interval=5/3, start=1.0, t_end=-1.0, zeromedian=true)],
    ),
    DatasetSpec(
        :FivePrime_gene,
        "trace",
        "data/5Prime_gene_enhancer/including_background/short",
        "gene";
        metadata=(sample=:FivePrime,),
        trace_specs=[(unit=1, interval=5/3, start=1.0, t_end=-1.0, zeromedian=true)],
    ),
]
```

The rate-family map tells the compiler which base rate positions belong to each family:

```julia
rate_families = Dict(
    :G => collect(1:4),
    :initiation => [5],
    :R => collect(6:8),
    :decay => [9],
    :noise => collect(10:13),
)
```

The recursive sharing object passed to `fit` is a named tuple:

```julia
noise_priors_by_group = Dict(
    :ThreePrime => [-0.018007, 0.214861, 1.256466, 0.405046],
    :FivePrime => [-0.004966, 0.118037, 0.607264, 0.223634],
)

sharing = (
    kind=:recursive,
    levels=grouping,
    datasets=datasets,
    transition_families=[:G, :initiation, :R, :decay],
    emission_families=[:noise],
    rate_families=rate_families,
    initial_values=Dict(:noise => noise_priors_by_group),
)
```

Then call `fit` as usual. The keyword is still `hierarchical` because that is the existing `fit` entry point for structured parameter layouts:

```julia
fit(;
    key="recursive-3Prime-5Prime-gene-3301",
    nchains=2,
    datatype="trace",
    datapath="",
    root="/path/to/project",
    resultfolder="/path/to/results",
    gene="MYC",
    cell="HBEC",
    datacond="gene",
    label="recursive-3Prime-5Prime-gene-3301",
    G=3,
    R=3,
    S=0,
    insertstep=1,
    transitions=([1, 2], [2, 1], [2, 3], [3, 2]),
    fittedparam=[collect(1:8); collect(10:13)],
    fixedeffects=tuple(),
    coupling=tuple(),
    grid=nothing,
    priormean=Float64[],
    priorcv=Float64[],
    noisepriors=noise_priors_by_group[:ThreePrime],
    decayrate=-1.0,
    probfn=prob_Gaussian,
    hierarchical=sharing,
    ratetype="ml",
    propcv_rate=0.05,
    propcv_noise=0.001,
    maxtime=100.0,
)
```

With two groups in this layout, the expanded base parameter count is:

```text
top:   G(4) + initiation(1) + decay(1) = 6
group: R(3) + noise(4) = 7 per group
total: 6 + 2*7 = 20
```

## Individual Noise Level

To fit noise separately for each trace, add an individual level and move `:noise` to it. This still does not add hyperparameters; it only changes how many noise copies are fit.

```julia
grouping_h = HierarchySpec(
    HierarchyNode(:top, children=[
        HierarchyNode(:group, key=:sample, levels=[:ThreePrime, :FivePrime], children=[
            HierarchyNode(:individual, key=:trace_id)
        ])
    ]);
    parameter_scope=Dict(
        :G => :top,
        :initiation => :top,
        :R => :group,
        :decay => :top,
        :noise => :individual,
    ),
)
```

The loader creates `trace_id` metadata when it expands folder-level `DatasetSpec` rows into per-trace assignments. In batch script names, use `-h` for this mode. No `-h` means no individual noise level.

## Multiple Levels

The grouping tree can have any number of nested levels. For example, one can share `G` globally, fit `R` by sample within perturbation, and fit noise by replicate:

```text
top
  perturbation
    sample
      replicate
```

```julia
grouping = HierarchySpec(
    HierarchyNode(:top, children=[
        HierarchyNode(:perturbation, key=:perturbation, levels=[:UNT, :KRAB, :VPR], children=[
            HierarchyNode(:sample, key=:sample, levels=[:ThreePrime, :FivePrime], children=[
                HierarchyNode(:replicate, key=:replicate)
            ])
        ])
    ]);
    parameter_scope=Dict(
        :G => :top,
        :initiation => :top,
        :R => :sample,
        :decay => :top,
        :noise => :replicate,
    ),
)
```

Rows must include metadata for every keyed level:

```julia
DatasetSpec(
    :UNT_ThreePrime_rep1_gene,
    "trace",
    "CRISPR/5Prime/UNT/rep1/short",
    "gene";
    metadata=(
        perturbation=:UNT,
        sample=:ThreePrime,
        replicate=:rep1,
    ),
)
```

The ownership key for a parameter is the path prefix down to that level. If `:R => :sample`, then the key includes both the parent perturbation and the sample:

```julia
(top=:top, perturbation=:UNT, sample=:ThreePrime)
```

This prevents `:ThreePrime` under `:UNT` from being accidentally merged with `:ThreePrime` under `:KRAB`.

## Cache Behavior

The recursive likelihood builds only the unique objects it needs:

- `transition_families` determines how many transition matrices `a` are cached.
- `emission_families` determines how many emission parameter groups are cached.
- Each trace is assigned to one transition group and one emission group.

For the no-`h` 3Prime/5Prime example, there are two transition groups and two emission groups. For the `-h` individual-noise case, there are still two transition groups, but emission groups are trace-specific.

## Script Generator Convention

Project script generators can expose this through a simple switch:

```julia
individual_noise=false  # no -h, noise is group-level
individual_noise=true   # -h, noise is individual trace-level
```

The generated fit script should print the chosen mode before fitting, for example:

```text
individual noise level: false
```

If the printed mode or filename does not match the intended setting, regenerate the script from a fresh Julia session before running the fit.

## What This Is Not

This interface does not currently create learned hyperparameters. For example, if `:R => :group`, then each group has its own fitted `R` copy:

```text
R_ThreePrime ~ prior
R_FivePrime  ~ prior
```

It does not fit:

```text
R_hyper_mean, R_hyper_cv
R_ThreePrime ~ distribution(R_hyper_mean, R_hyper_cv)
R_FivePrime  ~ distribution(R_hyper_mean, R_hyper_cv)
```

That distinction matters. This is a Bayesian grouped-parameter model with exact sharing where requested, not a hierarchical shrinkage model.
