# This file is part of StochasticGene.jl
#
# test.jl
#
# Two roles:
# (1) **Smoke tests** used by `test/runtests.jl`: implement scenarios as `test_*` functions
#     here so `runtests.jl` only calls them and applies `@test`; you can rerun any check
#     interactively (e.g. `test_fit_simrna()`). Supplementary suites may live under `test/*.jl`.
# (2) **Developer diagnostics** (e.g. detailed RG vs CoupledFull comparisons,
#     long correlation/simulation checks). These are NOT called from `runtests.jl`
#     and are intended for algorithm development and debugging.

# ════════════════════════════════════════════════════════════════════════════════════
# SECTION 1: UTILITY & HELPER FUNCTIONS
# ════════════════════════════════════════════════════════════════════════════════════
# Low-level infrastructure for model testing, simulator validation, and matrix operations.
# These functions support the higher-level tests that follow.

"""
    test_steadystatemodel(model::AbstractGMmodel, nhist)

Compares the chemical master solution to a Gillespie simulation for the steady-state mRNA distribution.

# Arguments
- `model`: An instance of `AbstractGMmodel`.
- `nhist`: The number of histogram bins for the simulation.

# Description
This function compares the steady-state mRNA distribution obtained from the chemical master equation solution to that obtained from a Gillespie simulation. It uses the rates and number of states from the provided model.

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: The steady-state mRNA distributions from the chemical master solution and the Gillespie simulation.
"""
function test_steadystatemodel(model::AbstractGMmodel, nhist)
    G = model.G
    r = model.rates
    g1 = steady_state(r[1:2*G], G - 1, nhist, model.nalleles)
    g2 = simulatorGM(r[1:2*G], G - 1, nhist, model.nalleles)
    return g1, g2
end

"""
    test_model(data::RNAOnOffData, model::AbstractGRSMmodel)

Simulates the RNA on-off model with splicing and compares it to the provided data.

# Arguments
- `data`: An instance of `RNAOnOffData` containing the observed data.
- `model`: An instance of `AbstractGRSMmodel` representing the model to be tested.

# Description
This function simulates the RNA on-off model with splicing using the provided model parameters and compares the simulated results to the observed data. It uses the `telegraphsplice0` function to perform the simulation.

# Returns
- `Nothing`: The function performs the simulation and comparison but does not return a value.
"""
function test_model(data::RNAOnOffData, model::AbstractGRSMmodel)
    telegraphsplice0(data.bins, data.nRNA, model.G - 1, model.R, model.rates, 1000000000, 1e-5, model.nalleles)
end

"""
    simDT_convert(v)

Convert dwell-time specs from 4-column layout to simulator layout.

Used by coupled dwell-time test helpers where only ON/OFF histogram pairs are
needed by the simulator path.
"""
simDT_convert(v) = [[s[1], s[3]] for s in v]

"""
    test_DT(; r, transitions, G, R, S, insertstep, onstates, dttype, bins)

Smoke-test dwell-time prediction for a single (uncoupled) model.

Builds dwell-time reporter/components and returns `predictedarray(...)` output.
"""
function test_DT(r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)])
    reporter, components = make_reporter_components_DT(transitions, G, R, S, insertstep, "", onstates, dttype)
    predictedarray(r, components, bins, reporter, dttype)
end

##### Developer diagnostics below (not used by `runtests.jl`) #####
# Swarm / run-spec helpers: `biowulf.jl` (`makeswarmfiles`, …). Combined-rate CSV drivers: `io.jl` (`create_combined_files*`, …).

"""
    test_sim(; r, transitions, G, R, S, insertstep, nhist, nalleles, onstates, bins, total, tol, ejectnumber)

Run an uncoupled stochastic simulation with dwell-time configuration.

Returns raw `simulator(...)` output for quick parity checks against analytical paths.
"""
function test_sim(r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nhist=150, nalleles=2, onstates=[Int[], [2, 3]], bins=[collect(5/3:5/3:200), collect(0.1:0.1:20)], total=10000000, tol=1e-6, ejectnumber=1)
    simulator(r, transitions, G, R, S, insertstep, nhist=nhist, nalleles=nalleles, onstates=onstates, bins=bins, totalsteps=total, tol=tol, ejectnumber=ejectnumber)
end

"""
    get_3unit_model_params()

Returns a NamedTuple with default parameters for the 3-unit coupled model with hidden unit.

This is a reference model used in [`test_fit_tracejoint_3unit`](@ref) featuring:
- Unit 1 (observable): G=2, R=1 (1 reporter, observable)
- Unit 2 (observable): G=2, R=1 (1 reporter, observable)
- Unit 3 (hidden): G=3, R=0 (no reporters, hidden from observation)
- Coupling: Unit 3 inhibits both units 1 and 2

# Returns
NamedTuple with fields:
- `transitions::Tuple`: State transitions for each unit
- `G, R, S, insertstep`: Model structure parameters
- `coupling::Tuple`: Coupling structure (unit 3 → units 1,2)
- `r::Vector{Float64}`: Default rate parameters (28 elements)
- `observed_units::Vector{Int}`: Which units are observable ([1, 2])
- `noiseparams::Vector{Int}`: Noise parameters per unit ([4, 4, 0])
- `trial_time::Float64`: Default simulation time per trial (720.0 minutes)
- `ntrials::Int`: Default number of trials (100)
- `lags::Vector{Int}`: Default lags for correlation (0:60)

# Example
```julia
params = get_3unit_model_params()
result = test_simulate_trials(ntrials=params.ntrials, trial_time=params.trial_time, lags=params.lags)
```
"""
function get_3unit_model_params()
    transitions = (([1, 2], [2, 1]), ([1, 2], [2, 1]), ([1, 2], [2, 1], [2, 3], [3, 2], [1, 3], [3, 1]))
    G = (2, 2, 3)
    R = (1, 1, 0)
    S = (0, 0, 0)
    insertstep = (1, 1, 0)
    coupling = ((1, 2, 3), [(3, 1, 1, 2), (3, 2, 2, 1)], [:inhibit, :inhibit])
    r_u1 = [0.1, 0.2, 0.5, 0.3, 1.0, 5, 3, 10, 2]
    r_u2 = [0.15, 0.25, 0.6, 0.4, 1.0, 5, 3, 10, 2]
    r_u3 = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3, 0.01]
    r_coupling = [-0.5, -0.7]
    r = vcat(r_u1, r_u2, r_u3, r_coupling)
    observed_units = [1, 2]
    noiseparams = [4, 4, 0]
    trial_time = 720.0
    ntrials = 1000
    lags = collect(0:60)
    return (
        transitions=transitions,
        G=G, R=R, S=S, insertstep=insertstep,
        coupling=coupling,
        r=r,
        observed_units=observed_units,
        noiseparams=noiseparams,
        trial_time=trial_time,
        ntrials=ntrials,
        lags=lags
    )
end

"""
    test_simulate_trials(; ntrials=100, trial_time=720.0, lags=collect(0:60), kwargs...)

Run [`simulate_trials`](@ref) with the default 3-unit coupled parameters from [`get_3unit_model_params`](@ref)
(two observed units, one hidden inhibitory unit, coupling, and per-unit `noiseparams`).

`interval=1.0` is passed through for frame timing and default `trace_specs`. Remaining `kwargs...` are forwarded to `simulate_trials`
(e.g. `warmupsteps`, `offset`, `correlation_algorithm`, `zeromedian`, `trace_specs`).

# Example
```julia
res = test_simulate_trials(ntrials=5, trial_time=80.0, warmupsteps=200)
```
"""
function test_simulate_trials(; ntrials=100, trial_time=720.0, lags=collect(0:60), kwargs...)
    params = get_3unit_model_params()
    simulate_trials(params.r, params.transitions, params.G, params.R, params.S,
                    params.insertstep, params.coupling, ntrials, trial_time, lags;
                    observed_units=params.observed_units, noiseparams=params.noiseparams,
                    interval=1.0, kwargs...)
end

"""
    test_cm(; r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates, dttype, bins, ejectnumber)

Compute CME predictions for dwell-time test configuration.

Builds `MComponents`/`MTComponents` and returns `predictedarray(...)`.
"""
function test_cm(r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nRNA=150, nalleles=2, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)], ejectnumber=1)
    reporter, tcomponents = make_reporter_components_DT(transitions, G, R, S, insertstep, "", onstates, dttype)
    mcomponents = MComponents(transitions, G, R, nRNA, r[num_rates(transitions, R, S, insertstep)], "", ejectnumber)
    components = MTComponents{typeof(mcomponents),typeof(tcomponents)}(mcomponents, tcomponents)
    predictedarray(r, components, bins, reporter, dttype, nalleles, nRNA)
end

"""
    test_CDT(; r, transitions, G, R, S, insertstep, onstates, dttype, bins, coupling)

Coupled dwell-time CME prediction helper (legacy RG stack defaults).

Returns `predictedarray(...)` using coupled reporter/components construction.
"""
function test_CDT(r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 1.0, 0.045, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 1.0, -0.5], transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1), onstates=[[Int[], Int[], [3], [3]], [Int[], Int[], [3], [3]]], dttype=[["ON", "OFF", "ONG", "OFFG"], ["ON", "OFF", "ONG", "OFFG"]], bins=[[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)], [collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]], coupling=((1, 2), [(1, 2, 2, 3)]))
    reporter, components = make_reporter_components_DT(transitions, G, R, S, insertstep, "", onstates, dttype, coupling)
    couplingindices = coupling_indices(transitions, R, S, insertstep, [0, 0], coupling, nothing)
    rates, couplingStrength = prepare_rates_coupled(r, [num_rates(transitions[1], R[1], S[1], insertstep[1]), num_rates(transitions[2], R[2], S[2], insertstep[2])], couplingindices)
    predictedarray(rates, couplingStrength, components, bins, reporter, dttype)
end

"""
    test_CDT_full(r, transitions, G, R, S, insertstep, dwell_specs, coupling)

Coupled dwell-time CME prediction. dwell_specs lists only observed units (hidden units omitted).
Builds full onstates/dttype/bins (placeholder for hidden units), runs predictedarray, returns histograms
only for units in dwell_specs. Coupled stack selection is controlled by `coupled_stack`.
"""
function test_CDT_full(r, transitions, G, R, S, insertstep, dwell_specs, coupling; splicetype="", coupled_stack::Symbol=:full)
    unit_model = coupling[1]
    n_units = length(unit_model)
    observed = StochasticGene.observed_units_from_dwell_specs(dwell_specs)
    placeholder = dwell_specs[1]
    onstates_full = Vector{typeof(placeholder.onstates)}(undef, n_units)
    dttype_full = Vector{typeof(placeholder.dttype)}(undef, n_units)
    bins_full = Vector{typeof(placeholder.bins)}(undef, n_units)
    for k in 1:n_units
        j = findfirst(s -> s.unit == unit_model[k], dwell_specs)
        if j !== nothing
            onstates_full[k] = dwell_specs[j].onstates
            dttype_full[k] = dwell_specs[j].dttype
            bins_full[k] = dwell_specs[j].bins
        else
            # Hidden unit: valid placeholder. R=0 requires explicit G-state onstates.
            dttype_full[k] = placeholder.dttype
            bins_full[k] = placeholder.bins
            onstates_full[k] = R[k] == 0 ? [[G[k]] for _ in 1:length(placeholder.dttype)] : placeholder.onstates
        end
    end
    if coupled_stack === :full
        reporter, components = make_reporter_components_DT(transitions, G, R, S, insertstep, splicetype, onstates_full, dttype_full, coupling)
    elseif coupled_stack === :legacy
        sojourn = sojourn_states(onstates_full, G, R, S, insertstep, dttype_full)
        components = TDCoupledComponents(coupling, transitions, G, R, S, insertstep, sojourn, dttype_full, splicetype)
        sojourn_c = coupled_states(sojourn, coupling, components, G)
        nonzeros = coupled_states(nonzero_rows(components), coupling, components, G)
        reporter = (sojourn_c, nonzeros)
    else
        throw(ArgumentError("coupled_stack must be :full or :legacy (got $(coupled_stack))"))
    end
    nrates = [num_rates(transitions[i], R[i], S[i], insertstep[i]) for i in eachindex(R)]
    couplingindices = coupling_indices(transitions, R, S, insertstep, zeros(Int, length(R)), coupling, nothing)
    rates, couplingStrength = prepare_rates_coupled(r, nrates, couplingindices)
    hists = predictedarray(rates, couplingStrength, components, bins_full, reporter, dttype_full)
    [hists[k] for k in 1:n_units if unit_model[k] in observed]
end

"""
    test_nullspace_solvers(M; nrep=10, verbose=true)

Benchmark and compare the QR-based and SVD-based nullspace solvers on
the given matrix `M` (typically a transition rate matrix). Returns a
NamedTuple with timings and basic accuracy diagnostics.
"""
function test_nullspace_solvers(M; nrep=10, verbose=true)
    # Warmup
    p_qr = StochasticGene.normalized_nullspace_qr(M)
    p_svd = StochasticGene.normalized_nullspace_svd(M)

    # Timings
    t_qr = 0.0
    t_svd = 0.0
    for _ in 1:nrep
        t_qr += @elapsed StochasticGene.normalized_nullspace_qr(M)
        t_svd += @elapsed StochasticGene.normalized_nullspace_svd(M)
    end
    t_qr /= nrep
    t_svd /= nrep

    # Accuracy diagnostics
    p_qr = StochasticGene.normalized_nullspace_qr(M)
    p_svd = StochasticGene.normalized_nullspace_svd(M)
    diff_l1 = sum(abs.(p_qr .- p_svd))
    diff_linf = maximum(abs.(p_qr .- p_svd))

    result = (t_qr=t_qr,
              t_svd=t_svd,
              speedup = t_svd > 0 ? t_qr / t_svd : Inf,
              diff_l1=diff_l1,
              diff_linf=diff_linf,
              sum_qr=sum(p_qr),
              sum_svd=sum(p_svd),
              p_qr=p_qr,
              p_svd=p_svd)

    if verbose
        println("QR nullspace:   time ≈ $(t_qr) s, sum=$(result.sum_qr)")
        println("SVD nullspace:  time ≈ $(t_svd) s, sum=$(result.sum_svd)")
        println("Relative speed (qr/svd): $(result.speedup)")
        println("‖qr - svd‖₁ = $(result.diff_l1), ‖qr - svd‖∞ = $(result.diff_linf)")
    end

    return result
end

"""
    test_set(; r, transitions, G, R, S, insertstep, coupling)

Return a compact default parameter bundle for simple 3-unit coupled tests.
"""
function test_set(;
    r=[0.38, 0.1, 0.23, 0.2, 1.0,
       0.45, 0.2, 0.43, 0.3, 1.0,
       0.11, 0.12, 0.13, 1.0, 
       -.9, -.7],
    transitions=(([1, 2], [2, 1]),
                 ([1, 2], [2, 1]),
                 ([1, 2], [2, 1])),
    G=(2, 2, 2), R=(1, 1, 0), S=(0, 0, 0), insertstep=(1, 1, 0),
    coupling=((1, 2, 3), [(3, 1, 1, 2), (3, 2, 2, 1)], [:free, :free]),
)
return r, transitions, G, R, S, insertstep, coupling
end


"""
    test_3unit(; r, transitions, G, R, S, insertstep, coupling, verbose=true)

Compare coupled transition matrices built by RG stack and Full stack for a 3-unit model.

Returns shape/magnitude agreement diagnostics and both dense matrices.
"""
function test_3unit(;
    r=[0.37, 0.1, 0.23, 0.21, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 1.0,
       0.71, 0.23, 0.43, 0.31, 0.2, 0.52, 0.31, 0.3, 0.86, 0.5, 1.0,
       0.11, 0.12, 0.13, 0.14, 0.15, 1.0, 
        -.9, -.7],
    transitions=(([1, 2], [2, 1], [2, 3], [3, 2]),
                 ([1, 2], [2, 1], [2, 3], [3, 2]),
                 ([1, 2], [2, 1], [2, 3], [3, 2])),
    G=(3, 3, 3), R=(3, 3, 0), S=(2, 2, 0), insertstep=(1, 1, 0),
    coupling=((1, 2, 3), [(3, 1, 1, 2), (3, 2, 2, 1)], [:free, :free]),
    verbose::Bool=true,
)
    # RG and full components for the 3-unit coupled case (legacy argument structure)
    comp_rg   = TCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
    comp_full = TCoupledFullComponents(coupling, transitions, G, R, S, insertstep, "")

    nrates = [num_rates(transitions[i], R[i], S[i], insertstep[i]) for i in eachindex(R)]
    couplingindices = coupling_indices(transitions, R, S, insertstep, zeros(Int, length(R)), coupling, nothing)
    rates, couplingStrength = prepare_rates_coupled(r, nrates, couplingindices)
    T_RG   = make_mat_TC(comp_rg,   rates, couplingStrength)

    couplingindices, targets = coupling_indices_full(transitions, R, S, insertstep, zeros(Int, length(R)), coupling, nothing)
    rates, coupling_rates = prepare_rates_coupled_full(r, nrates, couplingindices, targets)
    T_full = make_mat_TC(comp_full, rates, coupling_rates)


    # Matrix agreement diagnostics
    M_RG   = Matrix(T_RG)
    M_full = Matrix(T_full)
    same_shape = size(M_RG) == size(M_full)
    diff_max  = same_shape ? maximum(abs.(M_RG .- M_full)) : NaN
    matrices_match = same_shape && isapprox(M_RG, M_full; rtol=1e-12)

    if verbose
        println("3-unit T matrices: same_shape=$same_shape, max |T_RG - T_full| = $diff_max, match=$matrices_match")
    end

    return (same_shape = same_shape,
            diff_max = diff_max,
            matrices_match = matrices_match,
            M_RG=M_RG,
            M_full=M_full)
end

"""
    test_sim_coupled(; r, transitions, G, R, S, insertstep, onstates, dttype, bins, coupling, total, tol, verbose=false)

Run coupled simulator for dwell-time histograms with converted dwell specs.
"""
function test_sim_coupled(; r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 0.0, 0.2, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 0.0, -0.5], transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1), onstates=[[Int[], Int[], [2], [2]], [Int[], Int[], [2], [2]]], dttype=[["ON", "OFF", "ONG", "OFFG"], ["ON", "OFF", "ONG", "OFFG"]], bins=[[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)], [collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]], coupling=((1, 2), [(1, 2, 2, 2)]), total=10000000, tol=1e-6, verbose=false)
    simulator(r, transitions, G, R, S, insertstep, coupling=coupling, nhist=0, noiseparams=0, onstates=simDT_convert(onstates), bins=simDT_convert(bins), totalsteps=total, tol=tol, verbose=verbose)
end

"""
    sim_grid(; r, p, Ngrid, transitions, G, R, S, insertstep, totaltime, interval, ntrials)

Simulate grid-based traces (a_grid) for quick testing of grid HMM paths.
"""
function sim_grid(; r=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70], p=0.2, Ngrid=4, transitions=([1, 2], [2, 1]), G=2, R=1, S=1, insertstep=1, totaltime=1000.0, interval=1.0, ntrials=10)
    Nstate = num_rates(transitions, R, S, insertstep)
    a_grid = StochasticGene.make_a_grid(p, Ngrid)
    # traces = simulator(r, transitions, G, R, S, insertstep, traceinterval=interval, nhist=0, totaltime=totaltime, reporterfn=sum, a_grid=StochasticGene.make_a_grid(1.0, 4))[1]
    StochasticGene.simulate_trace_vector(r, transitions, G, R, S, insertstep, interval, totaltime, ntrials, a_grid=a_grid)
end

### functions used in runtest

"""
    test_compare(; r, transitions, G, R, S, insertstep, nRNA, nalleles, bins, total, tol, onstates, dttype)

Compare simulated and chemical master equation histograms for a given parameter set.

# Arguments
- `r`: Rate parameters.
- `transitions`, `G`, `R`, `S`, `insertstep`: Model structure.
- `nRNA`, `nalleles`: RNA and allele counts.
- `bins`: Histogram bins.
- `total`: Number of simulation steps.
- `tol`: Simulation tolerance.
- `onstates`, `dttype`: State and dwell time types.

# Returns
- Tuple of (chemical master histogram, array of simulated histograms).
"""

# ════════════════════════════════════════════════════════════════════════════════════
# SECTION 2: CME vs SIMULATOR COMPARISON TESTS
# ════════════════════════════════════════════════════════════════════════════════════
# Validate accuracy of chemical master equation against Gillespie simulator.


function test_compare(; r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nRNA=150, nalleles=2, bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)], total=10000000, tol=1e-6, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"], ejectnumber=1)
    hs = test_sim(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], total, tol, Int(ejectnumber))
    h = test_cm(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates, dttype, bins, ejectnumber)
    hs = StochasticGene.normalize_histogram.(hs)
    return h, make_array(hs)
end

"""
    test_fit_simrna_compare(; rtarget, transitions, G, nRNA, nalleles, fittedparam, fixedeffects, rinit, totalsteps, nchains, ejectnumber)

Generate synthetic RNA histogram data, fit the model, and compare predicted vs simulated histograms.

Returns `(hc, h)` where `hc` is the fitted model prediction and `h` is the raw simulator histogram.
"""
function test_fit_simrna_compare(; rtarget=[0.33, 0.19, 20.5, 1.0], transitions=([1, 2], [2, 1]), G=2, nRNA=100, nalleles=2, fittedparam=[1, 2, 3], fixedeffects=tuple(), rinit=[0.1, 0.1, 0.1, 1.0], totalsteps=100000, nchains=1, ejectnumber=1)
    h = simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=totalsteps, nalleles=nalleles, ejectnumber=ejectnumber)[1]
    data = RNAData{typeof(nRNA),typeof(h)}("", "", nRNA, h, 1.0,[])  # yield=1.0 (Float64, no inflation needed)
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0), fittedparam, fixedeffects, transitions, G, 0, 0, 0, "", nalleles, 10.0, Int[], rtarget[end], 0.02, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(1000000, 0, 20.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    hc = predictedfn(fits.parml, data, model)
    return hc, h
end

"""
    compare_data_model(data, model, parml)

Evaluate a model prediction against normalized RNA histogram data.

Returns `(predicted_hist, normalized_data_hist)`.
"""
function compare_data_model(data, model, parml)
    hc = predictedfn(parml, data, model)
    return hc, normalize_histogram(data.histRNA)
end

"""
    test_compare_coupling(; r, transitions, G, R, S, insertstep, onstates, dttype, bins, coupling, total, tol)

Compare simulated and chemical master equation histograms for coupled models.

# Arguments
- `r`: Rate parameters.
- `transitions`, `G`, `R`, `S`, `insertstep`: Model structure.
- `onstates`, `dttype`, `bins`: State and dwell time types, histogram bins.
- `coupling`: Coupling structure.
- `total`: Number of simulation steps.
- `tol`: Simulation tolerance.

# Returns
- Tuple of (chemical master histogram, array of simulated histograms).
"""
function test_compare_coupling(; r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 1.0, 0.45, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 1.0, 2.9], transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1), onstates=[[Int[], Int[], [3], [3]], [Int[], Int[], [3], [3]]], dttype=[["ON", "OFF", "ONG", "OFFG"], ["ON", "OFF", "ONG", "OFFG"]], bins=[[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)], [collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]], coupling=((1, 2), [(1, 3, 2, 4)], [:free]), total=1000000, tol=1e-6)
    hs = simulator(r, transitions, G, R, S, insertstep, coupling=coupling, nhist=0, noiseparams=0, onstates=simDT_convert(onstates), bins=simDT_convert(bins), totalsteps=total, tol=tol)
    h = test_CDT(r, transitions, G, R, S, insertstep, onstates, dttype, bins, coupling)
    for i in eachindex(hs)
        hs[i] = StochasticGene.normalize_histogram.(hs[i])
    end
    return make_array(vcat(h...)), make_array(vcat(hs...))
end

"""
    test_compare_coupling_full(; r, transitions, G, R, S, insertstep, onstates, dttype, bins, total, tol, verbose)

Compare RG-stack (`TCoupledComponents`) and Full-stack (`TCoupledFullComponents`) for reciprocal
coupling with mixed source types (connections in canonical α-first order):
- unit 2 → unit 1: R-step 2 specific (s = G₂+2), γ = r[21]
- unit 1 → unit 2: "any R occupied" sentinel (s = G₁+R₁+1), γ = r[22]

Checks:
1. T-matrix: RG vs Full (should match to machine precision). The "any R occupied" sentinel
   expands to all non-empty R-chain states via `classify_states`, so both stacks agree.
2. CME dwell-time histograms (Full-stack only) vs simulator. The RG dwell-time path uses
   a G×G source matrix and cannot represent R-state sentinels (s > G), so only the
   Full-stack (`TDCoupledFullComponents`) is used for histogram validation.

# Returns
NamedTuple with fields `T_match`, `T_diff`, `cme_full`, `sim`.
"""
function test_compare_coupling_full(;
    r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 1.0,
       0.45, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 1.0, 3.0, -0.7],
    transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])),
    G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1),
    onstates=[[Int[], Int[], [3], [3]], [Int[], Int[], [3], [3]]],
    dttype=[["ON", "OFF", "ONG", "OFFG"], ["ON", "OFF", "ONG", "OFFG"]],
    bins=[[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)],
          [collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]],
    total=1000000, tol=1e-6, verbose=true)
    # Canonical (α-first) order required so prepare_rates_sim and prepare_rates_coupled
    # assign γ values to the same connections:
    #   connection 1 (α=1, β=2): unit 2→unit 1, R-step 2 source (s = G₂+2), γ = r[21] = 3.0
    #   connection 2 (α=2, β=1): unit 1→unit 2, "any R occupied" sentinel (s = G₁+R₁+1), γ = r[22] = -0.7
    coupling = ((1, 2), [(2, G[2]+2, 1, 3), (1, G[1]+R[1]+1, 2, 3)], [:free, :free])

    nrates = [num_rates(transitions[i], R[i], S[i], insertstep[i]) for i in eachindex(R)]

    # 1. T-matrix: RG stack vs Full stack
    ci_rg = coupling_indices(transitions, R, S, insertstep, zeros(Int, length(R)), coupling, nothing)
    rates_rg, γ_rg = prepare_rates_coupled(r, nrates, ci_rg)
    ci_full, targets = coupling_indices_full(transitions, R, S, insertstep, zeros(Int, length(R)), coupling, nothing)
    rates_full, coupling_rates = prepare_rates_coupled_full(r, nrates, ci_full, targets)

    comp_rg   = TCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
    comp_full = TCoupledFullComponents(coupling, transitions, G, R, S, insertstep, "")
    M_RG   = Matrix(make_mat_TC(comp_rg,   rates_rg,   γ_rg))
    M_full = Matrix(make_mat_TC(comp_full, rates_full, coupling_rates))
    T_diff  = maximum(abs.(M_RG .- M_full))
    T_match = isapprox(M_RG, M_full; rtol=1e-12)
    verbose && println("T matrix: max|T_RG - T_full| = $T_diff, match = $T_match")

    # 2. CME dwell-time histograms (Full-stack only).
    # The RG dwell-time path uses a G×G source matrix (set_elements_Gs) that cannot
    # represent R-state source sentinels (s > G); only the Full-stack handles this.
    dwell_specs = [
        (unit=1, onstates=onstates[1], dttype=dttype[1], bins=bins[1]),
        (unit=2, onstates=onstates[2], dttype=dttype[2], bins=bins[2]),
    ]
    h_full = test_CDT_full(r, transitions, G, R, S, insertstep, dwell_specs, coupling; splicetype="")
    cme_full = make_array(vcat(h_full...))
    verbose && println("CME Full: $(length(cme_full)) histogram bins computed")

    # 3. CME Full vs simulator
    # warmupsteps ensures steady-state initial conditions match the CME's pss-based init
    hs = simulator(r, transitions, G, R, S, insertstep, coupling=coupling, nhist=0, noiseparams=0,
                   onstates=simDT_convert(onstates), bins=simDT_convert(bins), totalsteps=total, tol=tol,
                   warmupsteps=100000)
    for i in eachindex(hs)
        hs[i] = StochasticGene.normalize_histogram.(hs[i])
    end
    sim_vec = make_array(vcat(hs...))
    verbose && println("CME Full vs Simulator: len_cme=$(length(cme_full)), len_sim=$(length(sim_vec))")

    if verbose
        dtnames = ["ON", "OFF", "ONG", "OFFG"]
        for u in 1:length(h_full)
            for k in 1:length(h_full[u])
                cme_hist = h_full[u][k]
                sim_hist = hs[u][k]
                nm = k <= length(dtnames) ? dtnames[k] : "dt$k"
                bins_k = collect(1:length(cme_hist))
                mean_cme = sum(bins_k .* cme_hist)
                mean_sim = sum(bins_k .* sim_hist)
                println("Unit $u $nm: mean CME=$(round(mean_cme;digits=3)) sim=$(round(mean_sim;digits=3)) (ratio=$(round(mean_sim/mean_cme;digits=3)))")
                println("  CME head: ", round.(cme_hist[1:min(3,end)]; digits=5), " ... tail: ", round.(cme_hist[max(1,end-1):end]; digits=5))
                println("  sim head: ", round.(sim_hist[1:min(3,end)]; digits=5), " ... tail: ", round.(sim_hist[max(1,end-1):end]; digits=5))
                maxdiff = maximum(abs.(cme_hist .- sim_hist))
                println("  max|cme-sim|=$(round(maxdiff;digits=5))")
            end
        end
    end

    return (T_match=T_match, T_diff=T_diff,
            cme_full=cme_full, sim=sim_vec)
end

"""
    test_compare_3unit(; r, transitions, G, R, S, insertstep, total, tol, verbose=false)

Compare 3-unit coupled dwell-time CME predictions (full stack) against coupled simulator output.

The setup uses two observed units and one hidden unit with unit-specific coupling source states.
Returns `(cme_vec, sim_vec)` flattened arrays for direct numeric comparison.
"""
function test_compare_3unit(; r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 1.0, 0.45, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, -.9, -.7], transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1]), ([1, 2], [2, 1], [2, 3], [3, 2], [1, 3], [3, 1])), G=(2, 2, 3), R=(2, 1, 0), S=(1, 0, 0), insertstep=(1, 1, 0), total=1000000, tol=1e-6, verbose=false)
    coupling = ((1, 2, 3), [(3, 1, 1, 2), (3, 3, 2, 2)], [:free, :free])
    bins_u = [collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]
    # Unit 3 is hidden (no reporters, no onstates). Dwell specs only for observed units 1 and 2.
    dwell_specs = [
        (unit=1, onstates=[Int[], Int[], [2], [2]], dttype=["ON", "OFF", "ONG", "OFFG"], bins=bins_u),
        (unit=2, onstates=[Int[], Int[], [2], [2]], dttype=["ON", "OFF", "ONG", "OFFG"], bins=bins_u),
    ]
    hs = simulator_dwell_specs(r, transitions, G, R, S, insertstep, dwell_specs=dwell_specs, coupling=coupling, nhist=0, noiseparams=0, totalsteps=total, tol=tol)
    h = test_CDT_full(r, transitions, G, R, S, insertstep, dwell_specs, coupling; splicetype="")
    for i in eachindex(hs)
        hs[i] = StochasticGene.normalize_histogram.(hs[i])
    end
    # simulator_dwell_specs (coupled) returns one histogram per requested dt type (same shape as CME).
    cme_vec = make_array(vcat(h...))
    sim_vec = make_array(vcat(hs...))
    if verbose
        # Per-histogram comparison: label, length, sum, first 3 and last 2 bins (CME vs sim)
        dtnames = ["ON", "OFF", "ONG", "OFFG"]
        for u in 1:length(h)
            for k in 1:length(h[u])
                cme_hist = h[u][k]
                sim_hist = hs[u][k]
                nm = k <= length(dtnames) ? dtnames[k] : "dt$k"
                println("Unit $u $nm: len CME=$(length(cme_hist)) sim=$(length(sim_hist)) sum CME=$(round(sum(cme_hist); digits=6)) sim=$(round(sum(sim_hist); digits=6))")
                println("  CME head: ", round.(cme_hist[1:min(3,end)]; digits=5), " ... tail: ", round.(cme_hist[max(1,end-1):end]; digits=5))
                println("  sim head: ", round.(sim_hist[1:min(3,end)]; digits=5), " ... tail: ", round.(sim_hist[max(1,end-1):end]; digits=5))
            end
        end
    end
    return cme_vec, sim_vec
end

#= Feature/regression helpers moved to `test_features.jl`.

"""
    test_trace_full(; coupling, G, R, S, insertstep, transitions, rtarget, noisepriors,
                    interval, totaltime, ntrials, method)

Compare log-likelihoods from the legacy coupled stack (`TCoupledComponents`) and the
default coupled-full trace/HMM paths on the same simulated traces. The two paths must agree
to machine precision.

Returns a NamedTuple with fields `ll_rg`, `ll_full`, `diff`, `match`.
"""
function test_trace_full(;
    coupling=((1, 2), [(1, 2, 2, 1)]),
    G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1),
    transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])),
    rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 30, 100, 20,
             0.02, 0.05, 0.2, 0.2, 0.1, 50, 30, 100, 20, -0.5],
    noisepriors=([50, 30, 100, 20], [50, 30, 100, 20]),
    interval=1.0, totaltime=200.0, ntrials=3,
    method=Tsit5())

    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling,
                                  interval, totaltime, ntrials; verbose=false)
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], [0.0, 0.0], 1), Int[])
    rm = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, [5.0, 5.0], coupling)
    fittedparam = StochasticGene.set_fittedparam(Int[], data.label, transitions, R, S, insertstep, noisepriors, coupling, nothing)

    model_rg   = load_model(data, rtarget, rm, fittedparam, tuple(), transitions, G, R, S, insertstep, "",     1, 10.0, Int[], 1.0, 0.01, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing; coupled_stack=:legacy)
    model_full = load_model(data, rtarget, rm, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], 1.0, 0.01, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing; coupled_stack=:full)

    ll_rg   = loglikelihood(get_param(model_rg),   data, model_rg)[1]
    ll_full = loglikelihood(get_param(model_full), data, model_full)[1]
    diff = abs(ll_rg - ll_full)
    match = diff < 1e-8
    println("Trace ll: RG=$ll_rg, Full=$ll_full, |diff|=$diff, match=$match")
    return (ll_rg=ll_rg, ll_full=ll_full, diff=diff, match=match,
            components_rg=model_rg.components, components_full=model_full.components)
end

"""
    test_benchmark_trace_joint_fit_stacks(; nsamples, totaltime, ntrials, n_repeat, maxtime, ...)

Developer diagnostic: run [`run_mh`](@ref) on the **same** simulated coupled joint traces with
[`load_model`](@ref)`(...; coupled_stack=:legacy)` (`TCoupledComponents`, RG / Kronecker stack) vs
`:full` (`TCoupledFullComponents`). Prints mean wall time and allocations per stack (each stack runs
`n_repeat` independent fits). Not called from `runtests.jl`.

**Defaults** use **two units with `G=(3,3)`, `R=(3,3)`**, `S=(0,0)`, `insertstep=(1,1)`, and the usual
3-state promoter transitions **`([1,2],[2,1],[2,3],[3,2])`** per unit (must match `G`). Override
`G`, `R`, `transitions`, and **`rtarget`** (length `sum(num_rates per unit) + noise + ncoupling`) for
other layouts.

Use small `nsamples` / `maxtime` for a quick local check; increase for more stable timings.
[`test_trace_full`](@ref) checks likelihood parity; this helper compares **fit-loop** cost only.

# Returns
NamedTuple: `time_legacy`, `time_full`, `alloc_legacy`, `alloc_full`, `ratio_time_full_over_legacy`,
`fits_legacy`, `fits_full` (last MCMC result per stack).
"""
function test_benchmark_trace_joint_fit_stacks(;
    coupling=((1, 2), [(1, 2, 2, 1)], [:free]),
    G=(3, 3), R=(3, 3), S=(0, 0), insertstep=(1, 1),
    transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])),
    rtarget=Float64[
        0.03, 0.1, 0.5, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2,
        0.0, 0.1, 0.5, 0.15,
        0.03, 0.1, 0.5, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2,
        0.0, 0.1, 0.9, 0.2,
        -0.4,
    ],
    nsamples=1500,
    warmupsteps=0,
    annealsteps=0,
    totaltime=800.0,
    ntrials=4,
    fittedparam=Int[19],
    interval=1.0,
    noisepriors=([0., .1, 1., .1], [0., .1, 1., .1]),
    propcv=0.2,
    maxtime=120.0,
    method=Tsit5(),
    n_repeat=2,
    verbose=true,
    trace_specs=nothing,
)
    trace_specs_eff = trace_specs === nothing ? StochasticGene.default_trace_specs_for_coupled((interval, 1.0, -1.0), [false, false], 2) : trace_specs
    units = [spec.unit for spec in trace_specs_eff]
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling,
                                  interval, totaltime, ntrials; noiseparams=[4, 4])
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], fill(0.0, length(units)), 1), units)
    priormean = StochasticGene.set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors,
        StochasticGene.mean_elongationtime(rtarget, transitions, R), tuple(), coupling, nothing)
    rinit = rtarget
    options = StochasticGene.MHOptions(nsamples, warmupsteps, maxtime, 1.0)
    nr = StochasticGene.num_rates(transitions, R, S, insertstep)

    function bench_stack(stack::Symbol)
        time_acc = 0.0
        alloc_acc = 0
        fits_last = nothing
        for _ in 1:n_repeat
            model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep,
                "", 1, 10.0, Int[], rtarget[nr], propcv, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing;
                coupled_stack=stack)
            time_acc += @elapsed begin
                fits_last, _, _ = run_mh(data, model, options)
            end
            alloc_acc += @allocated run_mh(data, model, options)
        end
        return (time_acc / n_repeat, alloc_acc ÷ n_repeat, fits_last)
    end

    time_legacy, alloc_legacy, fits_legacy = bench_stack(:legacy)
    time_full, alloc_full, fits_full = bench_stack(:full)
    ratio = time_full / max(time_legacy, eps(Float64))
    if verbose
        println("test_benchmark_trace_joint_fit_stacks (nsamples=$nsamples, n_repeat=$n_repeat, ntrials=$ntrials):")
        println("  legacy (TCoupledComponents):  time=$(round(time_legacy; digits=4)) s, alloc=$alloc_legacy bytes")
        println("  full   (TCoupledFullComponents): time=$(round(time_full; digits=4)) s, alloc=$alloc_full bytes")
        println("  time_full / time_legacy = $(round(ratio; digits=3))")
    end
    return (time_legacy=time_legacy, time_full=time_full, alloc_legacy=alloc_legacy, alloc_full=alloc_full,
            ratio_time_full_over_legacy=ratio, fits_legacy=fits_legacy, fits_full=fits_full)
end

"""
    test_trace_specs(; coupling, G, R, S, insertstep, transitions, rtarget, noisepriors,
                     interval, totaltime, ntrials, method)

Test that `trace_specs` / `observed_units` limits the HMM likelihood to the observed
unit only. Uses the same 2-unit coupled setup as `test_trace_full` with `units=[1]`
(unit 2 is hidden). The key invariant: changing unit 2 noise params does NOT change
the log-likelihood when only unit 1 is observed via `data.units`.

Returns a NamedTuple with fields `ll`, `ll_perturbed`, `diff`, `match`.
"""
function test_trace_specs(;
    coupling=((1, 2), [(1, 2, 2, 1)]),
    G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1),
    transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])),
    rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 30, 100, 20,
             0.02, 0.05, 0.2, 0.2, 0.1, 50, 30, 100, 20, -0.5],
    noisepriors=([50, 30, 100, 20], [50, 30, 100, 20]),
    interval=1.0, totaltime=200.0, ntrials=3,
    method=Tsit5())

    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling,
                                  interval, totaltime, ntrials; verbose=false)

    # Only unit 1 is observed; unit 2 is hidden
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], [0.0, 0.0], 1), [1])
    rm = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, [5.0, 5.0], coupling)
    fittedparam = StochasticGene.set_fittedparam(Int[], data.label, transitions, R, S, insertstep, noisepriors, coupling, nothing)

    model = load_model(data, rtarget, rm, fittedparam, tuple(), transitions, G, R, S, insertstep,
                       "", 1, 10.0, Int[], 1.0, 0.01, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing)
    ll = loglikelihood(get_param(model), data, model)[1]

    # Locate unit 2 noise params and perturb them: LL must be invariant
    n1 = num_rates(transitions[1], R[1], S[1], insertstep[1])
    n2 = num_rates(transitions[2], R[2], S[2], insertstep[2])
    noise2_start = n1 + length(noisepriors[1]) + n2 + 1
    noise2_end   = noise2_start + length(noisepriors[2]) - 1
    rtarget2 = copy(rtarget)
    rtarget2[noise2_start:noise2_end] .*= 10.0

    model2 = load_model(data, rtarget2, rm, fittedparam, tuple(), transitions, G, R, S, insertstep,
                        "", 1, 10.0, Int[], 1.0, 0.01, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing)
    ll_perturbed = loglikelihood(get_param(model2), data, model2)[1]

    diff = abs(ll - ll_perturbed)
    match = diff < 1e-8
    println("trace_specs: ll=$ll, ll_perturbed_unit2=$ll_perturbed, |diff|=$diff, match=$match")
    return (ll=ll, ll_perturbed=ll_perturbed, diff=diff, match=match)
end

"""
    test_correlation_functions(; r, transitions, G, R, S, insertstep, coupling, lags, probfn, verbose)

Compute `correlation_functions` for the legacy (`coupled_stack=:legacy`) and default coupled-full
(`coupled_stack=:full`) paths and compare results. Because both stacks produce equivalent transition
matrices the outputs must agree to machine precision.

Returns a NamedTuple with fields `result_rg`, `result_full`, `diff_cc`, `diff_ac1`, `diff_ac2`,
`diff_ccON`, `match`.
"""
function test_correlation_functions(;
    coupling=((1, 2), [(1, 2, 2, 1)]),
    G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1),
    transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])),
    r=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 30, 100, 20,
       0.02, 0.05, 0.2, 0.2, 0.1, 50, 30, 100, 20, -0.5],
    lags=collect(0:1:20),
    probfn=prob_Gaussian,
    verbose=true)

    res_rg   = correlation_functions(r, transitions, G, R, S, insertstep, probfn, coupling, lags; splicetype="", coupled_stack=:legacy)
    res_full = correlation_functions(r, transitions, G, R, S, insertstep, probfn, coupling, lags; splicetype="", coupled_stack=:full)

    tau_rg, cc_rg, ac1_rg, ac2_rg, _, _, _, _, ccON_rg = res_rg[1:9]
    _,      cc_full, ac1_full, ac2_full, _, _, _, _, ccON_full = res_full[1:9]

    diff_cc   = maximum(abs.(cc_rg   .- cc_full))
    diff_ac1  = maximum(abs.(ac1_rg  .- ac1_full))
    diff_ac2  = maximum(abs.(ac2_rg  .- ac2_full))
    diff_ccON = maximum(abs.(ccON_rg .- ccON_full))
    match = diff_cc < 1e-8 && diff_ac1 < 1e-8 && diff_ac2 < 1e-8 && diff_ccON < 1e-8

    if verbose
        println("correlation_functions RG vs full:")
        println("  max|cc_rg   - cc_full|   = $diff_cc")
        println("  max|ac1_rg  - ac1_full|  = $diff_ac1")
        println("  max|ac2_rg  - ac2_full|  = $diff_ac2")
        println("  max|ccON_rg - ccON_full| = $diff_ccON")
        println("  match = $match")
    end
    return (result_rg=res_rg, result_full=res_full,
            diff_cc=diff_cc, diff_ac1=diff_ac1, diff_ac2=diff_ac2, diff_ccON=diff_ccON,
            match=match)
end

"""
    test_predict_traces(; rtarget, transitions, G, R, S, insertstep, coupling, noisepriors,
                         splicetype, interval, totaltime, ntrials, method)

Simulate coupled traces then compute predicted state sequences and observation distributions
via the full HMM decode path (`make_traces_dataframe`). Returns a NamedTuple with the
resulting DataFrames for both the RG-stack and the full-stack so the caller can compare.

No files are written; this is the "given the info, compute" equivalent of `write_traces_key`.
"""
function test_predict_traces(;
    coupling=((1, 2), [(1, 2, 2, 1)]),
    G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1),
    transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])),
    rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 30, 100, 20,
             0.02, 0.05, 0.2, 0.2, 0.1, 50, 30, 100, 20, -0.5],
    noisepriors=([50, 30, 100, 20], [50, 30, 100, 20]),
    interval=1.0, totaltime=200.0, ntrials=3,
    method=Tsit5())

    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling,
                                  interval, totaltime, ntrials; verbose=false)
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], [0.0, 0.0], 1), Int[])

    df_rg   = StochasticGene.make_traces_dataframe(data, rtarget, transitions, G, R, S, insertstep,
                                                    prob_Gaussian, 4, "",     true, false, coupling)
    df_full = StochasticGene.make_traces_dataframe(data, rtarget, transitions, G, R, S, insertstep,
                                                    prob_Gaussian, 4, "", true, false, coupling)

    println("test_predict_traces: df_rg rows=$(nrow(df_rg)), df_full rows=$(nrow(df_full))")
    return (df_rg=df_rg, df_full=df_full)
end

"""
    test_compare_RG_vs_Full(; r, dwell_specs, coupling, ntrials, rtol, verbose)

Run the same coupled dwell-time CME on the legacy coupled stack and the default coupled-full stack,
compare ON/OFF histograms and report time/memory. Uses 2-unit setup by default.
ONG/OFFG are not compared (legacy uses G-marginalized T; coupled-full uses expanded full T for all until that path exists).
Returns (h_RG, h_full, match::Bool, time_RG, time_full, alloc_RG, alloc_full).
"""
function test_compare_RG_vs_Full(;
    r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 1.0, 0.045, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 1.0, -0.5],
    transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])),
    G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1),
    dwell_specs=[(unit=1, onstates=[Int[], Int[], [3], [3]], dttype=["ON", "OFF", "ONG", "OFFG"], bins=[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]),
                 (unit=2, onstates=[Int[], Int[], [3], [3]], dttype=["ON", "OFF", "ONG", "OFFG"], bins=[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)])],
    coupling=((1, 2), [(1, 2, 2, 3)]),
    ntrials=2,
    rtol=1e-5,
    verbose=true)
    time_RG = 0.0
    time_full = 0.0
    alloc_RG = 0
    alloc_full = 0
    h_RG = nothing
    h_full = nothing
    for _ in 1:ntrials
        t = @elapsed h_RG = test_CDT_full(r, transitions, G, R, S, insertstep, dwell_specs, coupling; splicetype="", coupled_stack=:legacy)
        time_RG += t
        alloc_RG += @allocated test_CDT_full(r, transitions, G, R, S, insertstep, dwell_specs, coupling; splicetype="", coupled_stack=:legacy)
    end
    for _ in 1:ntrials
        t = @elapsed h_full = test_CDT_full(r, transitions, G, R, S, insertstep, dwell_specs, coupling; splicetype="", coupled_stack=:full)
        time_full += t
        alloc_full += @allocated test_CDT_full(r, transitions, G, R, S, insertstep, dwell_specs, coupling; splicetype="", coupled_stack=:full)
    end
    time_RG /= ntrials
    time_full /= ntrials
    alloc_RG = alloc_RG ÷ ntrials
    alloc_full = alloc_full ÷ ntrials
    # Compare ON/OFF only: full stack uses full T for all; RG uses G-marginalized T for ONG/OFFG
    h_RG_onoff = [[h_RG[α][1], h_RG[α][2]] for α in eachindex(h_RG)]
    h_full_onoff = [[h_full[α][1], h_full[α][2]] for α in eachindex(h_full)]
    cme_RG_onoff = make_array(vcat(h_RG_onoff...))
    cme_full_onoff = make_array(vcat(h_full_onoff...))
    match = isapprox(cme_RG_onoff, cme_full_onoff; rtol=rtol)
    if verbose
        println("RG   : time = $(round(time_RG; digits=6)) s, alloc = $alloc_RG bytes")
        println("Full : time = $(round(time_full; digits=6)) s, alloc = $alloc_full bytes")
        println("Match (rtol=$rtol): $match")
    end
    return h_RG, h_full, match, time_RG, time_full, alloc_RG, alloc_full
end

"""
    test_compare_T_matrix_RG_vs_Full(; r, coupling, rtol, verbose, uncoupled_only)

Build the full coupled transition matrix T using both the RG stack (Kronecker assembly)
and the Full stack (element expansion), then compare the two matrices. They should agree
to very high tolerance (default rtol=1e-12). Returns (T_RG, T_full, match::Bool).

If `uncoupled_only=true`, coupling_strength is zeroed so the comparison is uncoupled T only
(ensures expand_unit_elements_to_full matches Kronecker-sum construction).
"""
function test_compare_T_matrix_RG_vs_Full(;
    r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 1.0, 0.045, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 1.0, -0.5],
    transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])),
    G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1),
    coupling=((1, 2), [(1, 2, 2, 3)]),
    rtol=1e-12,
    verbose=true,
    uncoupled_only=false)
    nrates = [num_rates(transitions[i], R[i], S[i], insertstep[i]) for i in eachindex(R)]
    couplingindices, targets = coupling_indices_full(transitions, R, S, insertstep, zeros(Int, length(R)), coupling, nothing)
    rates, coupling_rates = prepare_rates_coupled_full(r, nrates, couplingindices, targets)
    if uncoupled_only
        coupling_rates = zeros(length(coupling_rates))
    end
    comp_rg   = TCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
    comp_full = TCoupledFullComponents(coupling, transitions, G, R, S, insertstep, "")
    # RG path uses raw coupling strengths (γ), full path uses precomputed coupling_rates (γ * base_rate)
    rg_coupling = uncoupled_only ? zeros(length(couplingindices)) : r[couplingindices]
    T_RG   = make_mat_TC(comp_rg,   rates, rg_coupling)
    T_full = make_mat_TC(comp_full, rates, coupling_rates)
    match = isapprox(Matrix(T_RG), Matrix(T_full); rtol=rtol)
    if verbose
        diff = maximum(abs.(Matrix(T_RG) - Matrix(T_full)))
        println("max |T_RG - T_full| = $diff" * (uncoupled_only ? " (uncoupled only)" : ""))
        println("Match (rtol=$rtol): $match")
    end
    return T_RG, T_full, match
end

"""
    test_compare_T_matrix_uncoupled_machine_precision(; kwargs...)

Run test_compare_T_matrix_RG_vs_Full with uncoupled_only=true and rtol=1e-14.
Assert that the uncoupled full matrix matches the RG (Kronecker-sum) matrix to machine precision.
"""
function test_compare_T_matrix_uncoupled_machine_precision(; kwargs...)
    T_RG, T_full, match = test_compare_T_matrix_RG_vs_Full(; uncoupled_only=true, rtol=1e-14, verbose=true, kwargs...)
    @assert match "Uncoupled T (full) must match T (RG) to machine precision; run test_compare_T_matrix_RG_vs_Full(uncoupled_only=true) for diagnostics."
    return T_RG, T_full
end

"""
    test_num_reporters_consistency(; G, R, insertstep, verbose)

Check that reporter count from CME (num_reporters_per_index) matches slice-based count
(num_reporters_from_r_slice) and simulator (num_reporters) for every state index z.

Builds the state matrix that corresponds to each z (R part = digit_vector(z, 2, R)) and
asserts num_reporters(state, ...) == num_reporters_per_index(z, ...). So we can tell:
- If slice vs per_index disagree → bug in transition_rate_functions.
- If state vs per_index disagree → bug in simulator's num_reporters or state encoding.

S=0 (base 2) only. Returns true if all checks pass.
"""
function test_num_reporters_consistency(; G=2, R=2, insertstep=1, verbose=false)
    R == 0 && return true
    S = 0
    base = 2
    nz = base^R
    for z in 1:nz
        d = StochasticGene.digit_vector(z, base, R)
        # Slice at and after insertstep (same as num_reporters_per_index)
        r_slice = d[insertstep:end]
        from_slice = StochasticGene.num_reporters_from_r_slice(r_slice)
        from_index = StochasticGene.num_reporters_per_index(z, R, insertstep, base, sum)
        if from_slice != from_index
            verbose && @warn "num_reporters mismatch: z=$z slice=$from_slice per_index=$from_index"
            return false
        end
        # State matrix encoding: R part = d (0/1). Simulator num_reporters must match.
        state = zeros(Int, G + R, 1)
        state[1, 1] = 1
        state[(G + 1):(G + R), 1] .= d
        from_sim = StochasticGene.num_reporters(state, 1, G, R, insertstep)
        if from_sim != from_index
            verbose && @warn "num_reporters sim vs CME: z=$z sim=$from_sim per_index=$from_index"
            return false
        end
    end
    verbose && @info "num_reporters consistent for G=$G R=$R insertstep=$insertstep ($nz states)"
    return true
end

"""
    diagnose_sim_vs_cme(; testfn=nothing, verbose=true)

Run a sim-vs-CME comparison and report where they disagree. Runs the comparison
(default: test_compare_3unit), prints per-histogram sums and a few bins so you can
see which histogram (ON/OFF/ONG/OFFG, which unit) is off, and returns (true, cme_vec, sim_vec).

Call as `diagnose_sim_vs_cme()` or `diagnose_sim_vs_cme(testfn=StochasticGene.test_compare, verbose=true)`.
"""
function diagnose_sim_vs_cme(; testfn=nothing, verbose=true)
    if testfn === nothing
        testfn = test_compare_3unit
    end
    verbose && println("=== Sim vs CME histograms ===")
    cme_vec, sim_vec = testfn(verbose=verbose)
    diff = abs.(cme_vec .- sim_vec)
    verbose && println("Max |CME - sim| = $(maximum(diff)); sum diff = $(sum(diff))")
    return true, cme_vec, sim_vec
end
=#

# Companion suites: feature/regression checks and developer benchmarks.
# Included here so these helpers remain available after `using StochasticGene`.
include("benchmarks.jl")
include("test_features.jl")

# ════════════════════════════════════════════════════════════════════════════════════
# SECTION 3: CORE INFERENCE TESTS
# ════════════════════════════════════════════════════════════════════════════════════
# Essential tests called by Pkg.test() (test/runtests.jl).
# These validate core package functionality: model fitting, AD gradients, and inference.
# Run in quick mode by default; extended tests available with FULL_TESTS=1 environment var.

"""
    test_fit_simrna(; rtarget, transitions, G, nRNA, nalleles, fittedparam, fixedeffects, rinit, totalsteps, nchains)

Fit a simulated RNA histogram using the provided parameters and compare to the target.

# Arguments
- `rtarget`: Target rate parameters.
- `transitions`, `G`, `nRNA`, `nalleles`: Model structure and counts.
- `fittedparam`, `fixedeffects`, `rinit`: Fitting options.
- `totalsteps`: Number of simulation steps.
- `nchains`: Number of MCMC chains.

# Returns
- Tuple of (fitted rates, target rates).
"""
function test_fit_simrna(; rtarget=[0.33, 0.19, 2.5, 1.0], transitions=([1, 2], [2, 1]), G=2, nRNA=100, nalleles=2, fittedparam=[1, 2, 3], fixedeffects=tuple(), rinit=[0.1, 0.1, 0.1, 1.0], totalsteps=100000, nchains=1)
    h = simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=totalsteps, nalleles=nalleles)[1]
    data = RNAData{typeof(nRNA),typeof(h)}("", "", nRNA, h, 1.0,[1])  # yield=1.0 (Float64, no inflation needed)
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0), fittedparam, fixedeffects, transitions, G, 0, 0, 0, "", nalleles, 10.0, Int[], rtarget[end], 0.02, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(1000000, 100000, 20.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    h = predictedfn(fits.parml, data, model)
    h, normalize_histogram(data.histRNA)
end

"""
    test_compare_yield(; r, transitions, G, R, S, insertstep, nRNA, nalleles, bins, total, tol, onstates, dttype, ejectnumber, yieldfactor)

Compare simulated and chemical master equation RNA histograms WITH yield factor (observation noise).

Like `test_compare` but includes yield < 1.0 for the RNA histogram only. Both the CME predictions 
AND simulated RNA data have observation loss applied via binomial sampling, then compared.
Dwell time histograms are returned unchanged.

# Arguments
- All arguments same as `test_compare`, plus:
- `yieldfactor::Float64`: Detection efficiency for RNA histogram (default 1.0). Values < 1.0 enable yield mode.

# Returns
- Tuple of (CME histograms [RNA_with_yield, dwell1, dwell2], simulated histograms [RNA_with_yield, dwell1, dwell2]).
"""
function test_compare_yield(; r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nRNA=150, nalleles=2, bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)], total=10000000, tol=1e-6, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"], ejectnumber=1, yieldfactor=0.5)
    hs = test_sim(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], total, tol, Int(ejectnumber))
    h = test_cm(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates, dttype, bins, ejectnumber)
    
    # Apply yield (technical loss) to RNA histogram only (first component)
    if yieldfactor < 1.0
        # Apply loss to simulated RNA histogram (hs is a tuple, so hs[1] is RNA)
        hs = (technical_loss(hs[1], yieldfactor, length(hs[1])), hs[2:end]...)
        
        # Apply loss to CME concatenated array - need to extract RNA part first
        # h is concatenated: [RNA..., ON..., OFF..., ONG...]
        h_rna_with_yield = technical_loss(h[1:nRNA], yieldfactor, nRNA)
        h = vcat(h_rna_with_yield, h[nRNA+1:end])
    end
    
    hs = StochasticGene.normalize_histogram.(hs)
    return h, make_array(hs)
end

"""
    test_fit_simrna_yield(; rtarget, transitions, G, nRNA, nalleles, fittedparam, fixedeffects, rinit, totalsteps, nchains, yieldfactor)

Fit synthetic RNA histogram WITH yield factor.

Generate synthetic RNA histogram data, apply yield-induced loss to create observed data, fit the model 
with that yield parameter, and compare predictions at both the true (lossless) and observed (lossy) levels.

# Arguments
- All arguments same as `test_fit_simrna`, plus:
- `yieldfactor::Float64`: Detection efficiency (default 0.5). Values < 1.0 enable yield mode
  where the model uses internal histogram expansion and loss matrix during fitting.

# Returns
- Named tuple with:
  - `full_pred`: Predicted full (lossless) histogram
  - `full_actual`: Actual full (true) simulated histogram (no loss applied)
  - `lossy_pred`: Predicted histogram after applying yieldfactor loss
  - `lossy_actual`: Observed histogram after applying yieldfactor loss to true sim data
"""
function test_fit_simrna_yield(; maxtime = 60.0, rtarget=[0.33, 0.19, 2.5, .1], transitions=([1, 2], [2, 1]), G=2, nRNA=100, nalleles=2, fittedparam=[1, 2, 3], fixedeffects=tuple(), rinit=[0.1, 0.1, 0.1, 1.0], totalsteps=1000000, nchains=1, yieldfactor=0.5)
    # Simulate without loss
    h_true = simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=totalsteps, nalleles=nalleles)[1]
    
    # Apply yield loss to create observed data
    if yieldfactor < 1.0
        h_observed_full = technical_loss(h_true, yieldfactor, nRNA)
        # Find the observed histogram size: where we reach 99% cumulative probability
        cdf_obs = cumsum(h_observed_full)
        nRNA_observed = findfirst(x -> x > 0.999, cdf_obs)
        h_observed = h_observed_full[1:nRNA_observed]
        nRNA_true = nhist_loss(nRNA_observed, yieldfactor)
        yield_tuple = (yieldfactor, nRNA_true)
        println("Yield=$yieldfactor: nRNA_true_original=$nRNA, nRNA_observed=$nRNA_observed, nRNA_true_estimated=$nRNA_true (expansion: $(nRNA_true/nRNA_observed)x)")
    else
        h_observed = h_true
        nRNA_observed = nRNA
        yield_tuple = yieldfactor
        println("No yield (yieldfactor=$yieldfactor)")
    end
    
    # Fit using the observed (lossy) data with yield parameter
    data = RNAData{typeof(nRNA_observed),typeof(h_observed)}("", "", nRNA_observed, h_observed, yield_tuple, [1])
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0), fittedparam, fixedeffects, transitions, G, 0, 0, 0, "", nalleles, 10.0, Int[], rtarget[end], 0.02, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    nRNA_true_actual = StochasticGene.get_nRNA_true(data.yield, data.nRNA)
    println("Model will use nRNA_true=$nRNA_true_actual for loss matrix")
    options = StochasticGene.MHOptions(5000000, 500000, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    
    # Get predictions at both levels
    h_pred_full = predictedfn(fits.parml, data, model)  # Prediction (normalized, loss matrix applied if yield < 1.0)
    
    # Compare using same scalings as loglikelihood does
    h_pred_lossy = h_pred_full  # Predictions from model
    h_actual_lossy = StochasticGene.datapdf(data)  # Data histogram normalized by datapdf
    
    return (
        full_actual=normalize_histogram(h_true),
        lossy_pred=h_pred_lossy,
        lossy_actual=h_actual_lossy
    )
end


"""
    test_fit_rna(; gene, G, nalleles, propcv, fittedparam, fixedeffects, transitions, rinit, datacond, datapath, label, root, nchains)

Fit a real RNA histogram using the provided parameters and compare to the data.

# Arguments
- `gene`, `G`, `nalleles`: Gene and model structure.
- `propcv`, `fittedparam`, `fixedeffects`, `transitions`, `rinit`: Fitting options.
- `datacond`, `datapath`, `label`, `root`: Data and file options.
- `nchains`: Number of MCMC chains.

# Returns
- Tuple of (predicted histogram, normalized data histogram).
"""
function test_fit_rna(; gene="CENPL", G=2, nalleles=2, propcv=0.05, fittedparam=[1, 2, 3], fixedeffects=tuple(), transitions=([1, 2], [2, 1]), rinit=[0.01, 0.1, 1.0, 0.01006327034802035], datacond="MOCK", datapath="data/HCT116_testdata", label="scRNA_test", root=".", nchains=1)
    data = load_data("rna", [], folder_path(datapath, root, "data"), label, gene, datacond, (), 1.0)
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, 1.0, [], 1.0), fittedparam, fixedeffects, transitions, 2, 0, 0, 1, "", nalleles, 10.0, Int[], rinit[end], propcv, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(100000, 100000, 60.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    h = predictedfn(fits.parml, data, model)
    h, normalize_histogram(data.histRNA)
end

"""
    test_fit_rnaonoff(; G, R, S, transitions, insertstep, rtarget, rinit, nsamples, nhist, nalleles, onstates, bins, fittedparam, propcv, priorcv, splicetype, nchains)

Fit a simulated RNA on-off histogram using the provided parameters and compare to the data.

# Arguments
- `G`, `R`, `S`, `transitions`, `insertstep`: Model structure.
- `rtarget`, `rinit`: Rate parameters.
- `nsamples`, `nhist`, `nalleles`, `onstates`, `bins`: Data and simulation options.
- `fittedparam`, `propcv`, `priorcv`, `splicetype`, `nchains`: Fitting options.

# Returns
- Tuple of (predicted histogram, data PDF).
"""
function test_fit_rnaonoff(; G=2, R=1, S=1, transitions=([1, 2], [2, 1]), insertstep=1, rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01], rinit=fill(0.01, num_rates(transitions, R, S, insertstep)), nsamples=100000, nhist=20, nalleles=2, onstates=Int[], bins=collect(1:1.0:200.0), fittedparam=collect(1:length(rtarget)-1), propcv=0.05, priorcv=10.0, splicetype="", nchains=1, yield=1.0)
    # OFF, ON, mhist = test_sim(rtarget, transitions, G, R, S, nhist, nalleles, onstates, bins)
    hs = simulator(rtarget, transitions, G, R, S, insertstep, nalleles=nalleles, nhist=nhist, onstates=onstates, bins=bins, totalsteps=10000000)
    hRNA = div.(hs[1], 30)
    data = StochasticGene.RNAOnOffData("test", "test", nhist, hRNA, bins, hs[2], hs[3], yield)
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], [], mean_elongationtime(rtarget, transitions, R)), fittedparam, tuple(), transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, rtarget[end], propcv, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 120.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    h = predictedfn(fits.parml, data, model)
    h, StochasticGene.datapdf(data)
end

"""
    test_fit_rnadwelltime(; rtarget, transitions, G, R, S, insertstep, nRNA, nsamples, nalleles, onstates, bins, dttype, fittedparam, propcv, priorcv, splicetype, maxtime, nchains)

Fit a simulated RNA dwell time histogram using the provided parameters and compare to the data.

# Arguments
- `rtarget`, `transitions`, `G`, `R`, `S`, `insertstep`: Model structure and rates.
- `nRNA`, `nsamples`, `nalleles`, `onstates`, `bins`, `dttype`: Data and simulation options.
- `fittedparam`, `propcv`, `priorcv`, `splicetype`, `maxtime`, `nchains`: Fitting options.

# Returns
- Tuple of (predicted histogram, data PDF).
"""
function test_fit_rnadwelltime(; rtarget=[0.038, 0.3, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.00231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nRNA=150, nsamples=100000, nalleles=2, onstates=[Int[], Int[], [2, 3], [2, 3]], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)], dttype=["ON", "OFF", "ONG", "OFFG"], fittedparam=collect(1:length(rtarget)-1), propcv=0.01, priorcv=10.0, splicetype="", maxtime=20.0, nchains=1, yield=1.0)
    h = test_sim(rtarget, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], 10000, 1e-2)
    data = RNADwellTimeData("test", "test", nRNA, h[1], bins, h[2:end], dttype, yield)
    rinit = StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], [], mean_elongationtime(rtarget, transitions, R))
    rinit=rtarget
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], [], mean_elongationtime(rtarget, transitions, R)), fittedparam, tuple(), transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, rtarget[end], propcv, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 1000, maxtime_seconds(maxtime), 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    h = predictedfn(fits.parml, data, model)
    println(fits.accept/fits.total)
    make_array(vcat(h...)), StochasticGene.datapdf(data)
    # lower = stats.qparam[1, :]
    # upper = stats.qparam[3, :]
    # return lower, rtarget[model.fittedparam], upper, fits, stats, measures
end

#= Feature/regression helpers now live in `test_features.jl`.

"""
    test_load_model_keyword_compatibility(; method=Tsit5())

Regression test for the `load_model(...; coupled_stack=...)` upstream path across non-trace
histogram data. The keyword is meaningful only for trace stacks, but RNA / ON-OFF / RNA+dwell
paths must accept the same upstream call without throwing a `MethodError`.
"""
function test_load_model_keyword_compatibility(; method=Tsit5())
    checks = Bool[]

    for stack in (:full, :legacy)
        transitions_rna = ([1, 2], [2, 1])
        r_rna = [0.02, 0.1, 0.5, 0.2]
        data_rna = RNAData("compat", "test", 10, ones(11), 1.0)
        rm_rna = prior_ratemean(transitions_rna, 1, 0, 1, r_rna[end], Float64[], mean_elongationtime(r_rna, transitions_rna, 1))
        model_rna = load_model(data_rna, r_rna, rm_rna, collect(1:length(r_rna)-1), tuple(), transitions_rna, 2, 1, 0, 1, "", 1, 10.0, Int[], r_rna[end], 0.05, prob_Gaussian, Float64[], method, tuple(), tuple(), nothing; coupled_stack=stack)
        push!(checks, model_rna.components isa MComponents)

        transitions_onoff = ([1, 2], [2, 1])
        r_onoff = [0.02, 0.1, 0.5, 0.2, 0.1, 0.01]
        data_onoff = RNAOnOffData("compat", "test", 10, ones(11), collect(0.1:0.1:2.0), ones(20), ones(20), 1.0)
        rm_onoff = prior_ratemean(transitions_onoff, 1, 1, 1, r_onoff[end], Float64[], mean_elongationtime(r_onoff, transitions_onoff, 1))
        model_onoff = load_model(data_onoff, r_onoff, rm_onoff, collect(1:length(r_onoff)-1), tuple(), transitions_onoff, 2, 1, 1, 1, "", 1, 10.0, Int[], r_onoff[end], 0.05, prob_Gaussian, Float64[], method, tuple(), tuple(), nothing; coupled_stack=stack)
        push!(checks, model_onoff.components isa MTAIComponents)

        transitions_dt = ([1, 2], [2, 1], [2, 3], [3, 2])
        r_dt = [0.038, 0.3, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.00231]
        onstates_dt = [Int[], Int[], [2, 3], [2, 3]]
        bins_dt = [collect(0.1:0.1:2.0) for _ in 1:4]
        dwell_dt = [ones(length(bins_dt[1])) for _ in 1:4]
        data_dt = RNADwellTimeData("compat", "test", 20, ones(21), bins_dt, dwell_dt, ["ON", "OFF", "ONG", "OFFG"], 1.0)
        rm_dt = prior_ratemean(transitions_dt, 2, 2, 1, r_dt[end], Float64[], mean_elongationtime(r_dt, transitions_dt, 2))
        model_dt = load_model(data_dt, r_dt, rm_dt, collect(1:length(r_dt)-1), tuple(), transitions_dt, 3, 2, 2, 1, "", 2, 10.0, onstates_dt, r_dt[end], 0.05, prob_Gaussian, Float64[], method, tuple(), tuple(), nothing; coupled_stack=stack)
        push!(checks, model_dt.components isa MTComponents)
    end

    all(checks)
end

struct MHRejectDummyModel <: AbstractGeneTransitionModel end
struct MHRejectDummyData <: AbstractExperimentalData end

num_rates(::MHRejectDummyModel) = 1
get_rates(param, ::MHRejectDummyModel, inverse::Bool=true) = param
logprior(param, ::MHRejectDummyModel) = 0.0
loglikelihood(param, ::MHRejectDummyData, ::MHRejectDummyModel) = throw(DimensionMismatch("synthetic invalid likelihood"))

"""
    test_mhstep_rejects_invalid_likelihood()

Regression test: if likelihood evaluation fails for a proposed point (e.g. malformed
RNA dwell-time solve output), Metropolis-Hastings should reject the proposal rather than
crash the run.
"""
function test_mhstep_rejects_invalid_likelihood()
    model = MHRejectDummyModel()
    data = MHRejectDummyData()
    param = [0.1]
    d = proposal_dist(param, 0.01, model)
    accept, logpredictions, param_out, ll_out, prior_out, _ = mhstep([0.0], param, 0.0, 0.0, d, 0.01, model, data, 1.0)
    accept == 0 && logpredictions == [0.0] && param_out == param && ll_out == 0.0 && prior_out == 0.0
end
=#

"""
    compute_quantile_bounds(qparam_row1, qparam_row3, n_accepts, n_total)

Expand quantile-based bounds to account for sampling variability from limited MCMC accepts.
With low acceptance rates, quantile estimates are uncertain and need wider margins.

# Arguments
- `qparam_row1`: Vector of 2.5% quantiles (lower bounds)
- `qparam_row3`: Vector of 97.5% quantiles (upper bounds)
- `n_accepts`: Number of accepted samples (from fits.accept)
- `n_total`: Total number of samples (from fits.total)

# Returns
- Tuple of (expanded_lower, expanded_upper) bounds

# Notes
Expansion factor scales with acceptance rate: lower acceptance → wider bounds.
Accounts for increased noise in quantile estimation with small effective sample sizes.
"""
function compute_quantile_bounds(qparam_row1, qparam_row3, n_accepts, n_total)
    # Effective sample size based on acceptance rate
    # Expansion factor increases as n_accepts decreases (more noisy quantiles)
    # At 283 accepts: ~13%, at 500: ~9%, at 1000: ~6%, at 2000: ~4%
    expansion_factor = 1.0 + 2.0 / sqrt(max(n_accepts, 1))
    
    lower = qparam_row1 ./ expansion_factor
    upper = qparam_row3 .* expansion_factor
    return lower, upper
end

"""
    test_fit_trace(; G, R, S, insertstep, transitions, rtarget, rinit, nsamples, onstates, totaltime, ntrials, fittedparam, propcv, cv, interval, noisepriors, nchains)

Fit a simulated trace dataset using the provided parameters and compare to the target.

# Arguments
- `G`, `R`, `S`, `insertstep`, `transitions`: Model structure.
- `rtarget`, `rinit`: Rate parameters.
- `nsamples`, `onstates`, `totaltime`, `ntrials`: Data and simulation options.
- `fittedparam`, `propcv`, `cv`, `interval`, `noisepriors`, `nchains`: Fitting options.

# Returns
- Tuple of (fitted rates, target rates).
"""
function test_fit_trace(; traceinfo=(1.0, 1.0, -1, 1.0, 0.5), G=2, R=2, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.01, 0.01, 0.1, 0.15, 0.1, 1.0, 50, 5, 50, 5], nsamples=10000, totaltime=4000.0, ntrials=10, fittedparam=[1:num_rates(transitions, R, S, insertstep)-1; num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+1], propcv=0.05, cv=100.0, noisepriors=[0.0, 0.1, 1.0, 0.1], nchains=1, zeromedian=true, maxtime=60.0, initprior=0.1)
    tracer = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, traceinfo[1], totaltime, ntrials)
    trace, tracescale = zero_median(tracer, zeromedian)
    nframes = round(Int, mean(size.(trace, 1)))  #mean number of frames of all traces
    if length(traceinfo) > 3 && traceinfo[4] != 1.0
        weight = set_trace_weight(traceinfo)
        background = set_trace_background(traceinfo)
    else
        weight = 0.0
        background = 0.0
    end
    data = TraceData{String,String,Tuple}("trace", "gene", traceinfo[1], (trace, background, weight, nframes, tracescale), Int[])
    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors, mean_elongationtime(rtarget, transitions, R), tuple(), tuple(), nothing)
    elongationtime = mean_elongationtime(rtarget, transitions, R)
    priormean = [fill(0.01, length(transitions)); initprior; fill(R / elongationtime, R); fill(0.05, max(0, S - insertstep + 1)); 1.0; noisepriors]
    priorcv = [fill(1.0, length(transitions)); 0.1; fill(0.1, R); fill(0.1, max(0, S - insertstep + 1)); 1.0; [0.5, 0.5, 0.1, 0.1]]
    rinit = isempty(tuple()) ? set_rinit(rtarget, priormean) : set_rinit(rtarget, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), tuple(), nothing)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, priorcv, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, Tsit5(), tuple(), tuple(), nothing, zeromedian)
    options = StochasticGene.MHOptions(nsamples, nsamples, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    println(fits.accept," ",fits.total)
    lower, upper = compute_quantile_bounds(stats.qparam[1, 1:4], stats.qparam[3, 1:4], fits.accept, fits.total)
    return lower, rtarget[1:4], upper
end

function test_fit_trace_compare(; traceinfo=(1.0, 1.0, -1, 1.0, 0.5), G=2, R=2, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.05, 0.2, 0.1, 0.15, 0.1, 1.0, 50, 5, 50, 5], nsamples=10000000, onstates=Int[], totaltime=1000.0, ntrials=100, fittedparam=[1:num_rates(transitions, R, S, insertstep)-1; num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+1], propcv=0.01, cv=100.0, noisepriors=[0.0, 0.2, 0.9, 0.1], nchains=1, zeromedian=true, maxtime=100.0, initprior=0.1)
    tracer = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, traceinfo[1], totaltime, ntrials)
    trace, tracescale = zero_median(tracer, zeromedian)
    nframes = round(Int, mean(size.(trace, 1)))  #mean number of frames of all traces
    if length(traceinfo) > 3 && traceinfo[4] != 1.0
        weight = set_trace_weight(traceinfo)
        background = set_trace_background(traceinfo)
    else
        weight = 0.0
        background = 0.0
    end
    data = TraceData{String,String,Tuple}("trace", "gene", traceinfo[1], (trace, background, weight, nframes, tracescale), Int[])
    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors, mean_elongationtime(rtarget, transitions, R), tuple(), tuple(), nothing)
    elongationtime = mean_elongationtime(rtarget, transitions, R)
    priormean = [fill(0.01, length(transitions)); initprior; fill(R / elongationtime, R); fill(0.05, max(0, S - insertstep + 1)); 1.0; noisepriors]
    priorcv = [fill(1.0, length(transitions)); 0.1; fill(0.1, R); fill(0.1, max(0, S - insertstep + 1)); 1.0; [0.5, 0.5, 0.1, 0.1]]
    rinit = isempty(tuple()) ? set_rinit([], priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), tuple(), nothing)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, priorcv, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, Tsit5(), tuple(), tuple(), nothing, zeromedian)
    options = StochasticGene.MHOptions(nsamples, 0, maxtime_seconds(maxtime), 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    return fits, stats, measures, data, model, options

end

"""
    test_fit_trace_hierarchical(; G, R, S, insertstep, transitions, rtarget, rinit, nsamples, onstates, totaltime, ntrials, fittedparam, propcv, cv, interval, noisepriors, hierarchical, method, maxtime, nchains)

Fit a simulated trace dataset using a hierarchical model and compare to the target.

# Arguments
- `G`, `R`, `S`, `insertstep`, `transitions`: Model structure.
- `rtarget`, `rinit`: Rate parameters.
- `nsamples`, `onstates`, `totaltime`, `ntrials`: Data and simulation options.
- `fittedparam`, `propcv`, `cv`, `interval`, `noisepriors`, `hierarchical`, `method`, `maxtime`, `nchains`: Fitting options.
- `hierarchical_sim_seed`: RNG seed for per-trial hierarchical rates `rh` (`nothing` uses the global RNG).
- `mh_seed`: `Random.seed!` before [`run_mh`](@ref) for reproducible chains (`nothing` skips).
- `warmupsteps`: Passed to [`MHOptions`](@ref) (with `nsamples` as `samplesteps`).

# Returns
- Tuple `(lower, target, upper)` of 2.5% / true / 97.5% quantiles for the first four rates.
"""
function test_fit_trace_hierarchical(; traceinfo=(1.0, 1.0, -1, 1.0, 0.5), G=2, R=2, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.01, 0.01, 0.1, 0.1, 0.1, 1.0, 5, 0.5, 20, 1], nsamples=10000, onstates=Int[], totaltime=4000.0, ntrials=200, fittedparam=[1, 2, 3, 4], propcv=0.01, noisepriors=[0.0, 0.1, 1.0, 0.1], hierarchical=(2, [7], tuple()), method=(Tsit5(), true), maxtime=120.0, nchains=1, zeromedian=true)
    rh = 5.0 .+ 0.5 * randn(ntrials)
    tracer = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, traceinfo[1], totaltime, ntrials, hierarchical=(6, rh))
    trace, tracescale = zero_median(tracer, zeromedian)
    nframes = round(Int, mean(size.(trace, 1)))  #mean number of frames of all traces
    if length(traceinfo) > 3 && traceinfo[4] != 1.0
        weight = set_trace_weight(traceinfo)
        background = set_trace_background(traceinfo)
    else
        weight = 0.0
        background = 0.0
    end
    data = TraceData{String,String,Tuple}("trace", "gene", traceinfo[1], (trace, background, weight, nframes, tracescale), Int[])
    # trace, tracescale = zero_median(tracer, zeromedian)
    # data = StochasticGene.TraceData("trace", "test", interval, (trace, [], 0.0, 1))(trace, background, weight, nframes, tracescale)
    rinit = []
    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors, mean_elongationtime(rtarget, transitions, R), hierarchical, tuple(), nothing)
    rinit = isempty(hierarchical) ? set_rinit(rinit, priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), tuple(), nothing; nhypersets=hierarchical[1])
    fittedparam = set_fittedparam(fittedparam, "trace", transitions, R, S, insertstep, noisepriors, tuple(), tuple())
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 0.1, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, hierarchical, tuple(), nothing, zeromedian)
    options = StochasticGene.MHOptions(nsamples, nsamples, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    println(fits.accept," ",fits.total)
    lower, upper = compute_quantile_bounds(stats.qparam[1, 1:4], stats.qparam[3, 1:4], fits.accept, fits.total)
    return lower, rtarget[1:4], upper
end

"""
    test_fit_tracejoint(; coupling, G, R, S, insertstep, transitions, rtarget, rinit, nsamples, onstates, totaltime, ntrials, fittedparam, propcv, cv, interval, noisepriors, maxtime, method)

Fit a simulated joint trace dataset for coupled models and compare to the target.

# Arguments
- `coupling`, `G`, `R`, `S`, `insertstep`, `transitions`: Model structure and coupling.
- `rtarget`, `rinit`: Rate parameters.
- `nsamples`, `onstates`, `totaltime`, `ntrials`: Data and simulation options.
- `fittedparam`, `propcv`, `cv`, `interval`, `noisepriors`, `maxtime`, `method`: Fitting options.

# Returns
- Tuple of (fitted rates, target rates).

# Notes
- Delegates to [`test_fit_tracejoint_full`](@ref): **full-stack** coupled traces (`TCoupledFullComponents`,
  coupled-full trace/HMM path, with `trace_specs` / `data.units` for observed units (same path as production fits).
"""
function test_fit_tracejoint(; kwargs...)
    test_fit_tracejoint_full(; kwargs...)
end


"""
    test_fit_tracejoint_full(; coupling, G, R, S, insertstep, transitions, rtarget, rinit,
                              nsamples, trace_specs, totaltime, ntrials, fittedparam,
                              propcv, interval, noisepriors, maxtime, method)

Coupled-full (`TCoupledFullComponents`) joint-trace test. Uses `trace_specs`
(default: [`default_trace_specs_for_coupled`](@ref)) so each spec has `unit`, `interval`, `start`,
`t_end`, `zeromedian`. Observed-unit indices are stored in `data.units` for correct HMM emission.

Fits the coupling strength (last parameter) and checks that `rtarget[fittedparam]` falls
within the posterior credible interval.

Returns `(lower, target, upper)` where `lower` and `upper` are the 2.5/97.5 percentiles.
"""
function test_fit_tracejoint_full(;
    coupling=((1, 2), [(1, 2, 2, 1)], [:free]),
    G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1),
    transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])),
    rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 1.0, 0., .05, 1.0, .05,
             0.03, 0.1, 0.5, 0.2, 1.0, 0.0, 0.05, 1.0, .05, -0.4],
    rinit=Float64[],
    nsamples=10000,
    trace_specs=nothing,
    totaltime=2000.0, ntrials=10,
    fittedparam=Int[21],
    propcv=0.2, interval=1.0,
    noisepriors=([0., .1, 1., .1], [0., .1, 1., .1]),
    maxtime="60m", method=Tsit5())

    trace_specs_eff = trace_specs === nothing ? StochasticGene.default_trace_specs_for_coupled((interval, 1.0, -1.0), [false, false], 2) : trace_specs
    units = [spec.unit for spec in trace_specs_eff]
    # Coupled simulator uses the legacy path; inference uses `TCoupledFullComponents`.
    # Per-unit noise counts: one entry per unit (hidden units with R=0 use 0).
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling,
                                  interval, totaltime, ntrials; noiseparams=[4, 4])
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], fill(0.0, length(units)), 1), units)
    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors,
                              mean_elongationtime(rtarget, transitions, R), tuple(), coupling, nothing)
    rinit = rtarget
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep,
                       "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)],
                       propcv, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing)
    options = StochasticGene.MHOptions(nsamples, 100, maxtime_seconds(maxtime), 1.0)
    fits, stats, measures = run_mh(data, model, options)
    println(fits.accept, " ", fits.total)
    lower, upper = compute_quantile_bounds(stats.qparam[1, :], stats.qparam[3, :], fits.accept, fits.total)
    return lower, rtarget[fittedparam], upper
end

#= Feature-level AD/NUTS/ADVI smoke tests moved to `test_features.jl`.

"""
    test_normalized_nullspace_augmented_pullback_fd(; n, ε, rtol, atol)

Finite-difference check for [`pullback_normalized_nullspace_augmented`](@ref). Intended for
`test/runtests.jl` and interactive runs. Returns `true` if all checks pass.
"""
function test_normalized_nullspace_augmented_pullback_fd(; n::Int=5, ε::Float64=1e-7, rtol::Float64=0.02, atol::Float64=1e-4)
    function example_row_sum_zero_sparse(n::Int)
        M = zeros(n, n)
        for i in 1:n, j in 1:n
            i == j && continue
            M[i, j] = 0.1 + 0.9 * mod1(i + 2j, 11) / 11
        end
        for i in 1:n
            M[i, i] = -sum(M[i, k] for k in 1:n if k != i)
        end
        sparse(M)
    end
    M = example_row_sum_zero_sparse(n)
    ȳ = randn(n)
    p, M̄ = pullback_normalized_nullspace_augmented(M, ȳ)
    length(p) == n || return false
    isapprox(sum(p), 1.0; atol=1e-9) || return false
    Ii, Jj, V = findnz(M)
    for k in eachindex(V)
        M2 = SparseMatrixCSC(copy(M))
        M2.nzval[k] += ε
        p2 = normalized_nullspace_augmented(M2)
        fd = (dot(ȳ, p2) - dot(ȳ, p)) / ε
        isapprox(M̄.nzval[k], fd; rtol=rtol, atol=atol) || return false
    end
    return true
end

"""
    test_trace_specs_utilities()

Checks [`n_observed_trace_units`](@ref) and [`default_trace_specs_for_coupled`](@ref) invariants
used by `test/runtests.jl`. Returns `true` if all checks pass.
"""
function test_trace_specs_utilities()
    c2 = ((1, 2), NTuple{4,Int}[], Symbol[])
    n_observed_trace_units(c2) == 2 || return false
    c3 = ((1, 2, 3), NTuple{4,Int}[], Symbol[])
    n_observed_trace_units(c3) == 2 || return false
    sp = default_trace_specs_for_coupled((1.0, 1.0, -1.0), [true, true], 2)
    return length(sp) == 2 && sp[1].unit == 1 && sp[2].unit == 2 && sp[1].interval == 1.0
end

function _test_twostate_generator(θ::AbstractVector{T}) where {T}
    length(θ) == 2 || throw(ArgumentError("expected length-2 rate vector"))
    θ1, θ2 = θ[1], θ[2]
    sparse(Int[1, 2, 1, 2], Int[1, 1, 2, 2], T[-θ2, θ2, θ1, -θ1], 2, 2)
end

function _test_central_grad(f, x::AbstractVector{<:Real}; ε::Float64=1e-6)
    g = Vector{Float64}(undef, length(x))
    δ = zeros(length(x))
    for i in eachindex(x)
        δ[i] = ε
        g[i] = (f(x .+ δ) - f(x .- δ)) / (2ε)
        δ[i] = 0.0
    end
    return g
end

"""
    test_ad_gradient_smoke(; ε=1e-6, rtol=0.02, atol=1e-3)

Lightweight AD installation check used by `test/runtests.jl`.
Verifies a few small augmented steady-state derivatives against central finite differences.
"""
function test_ad_gradient_smoke(; ε::Float64=1e-6, rtol::Float64=0.02, atol::Float64=1e-3)
    θ0 = Float64[1.2, 0.7]
    f(θ) = sum(normalized_nullspace_augmented(_test_twostate_generator(θ)))
    z = Zygote.gradient(f, θ0)[1]
    fd = _test_central_grad(f, θ0; ε=ε)
    isapprox(z, fd; rtol=rtol, atol=atol) || return false

    w = Float64[0.2, 0.8]
    g(θ) = dot(w, steady_state_vector(_test_twostate_generator(θ); solver=:augmented))
    zg = Zygote.gradient(g, θ0)[1]
    fdg = _test_central_grad(g, θ0; ε=ε)
    isapprox(zg, fdg; rtol=rtol, atol=atol) || return false

    return true
end

"""
    test_get_rates_ad_consistency(; ε=1e-4)

Smoke-test the AD-friendly rate reconstruction helpers used by histogram and trace likelihoods.
Returns `true` when fixed effects, `get_rates_ad`, and `predictedfn(...; rates_fn=get_rates_ad)`
match their non-AD counterparts on short GM / GRSM examples.
"""
function test_get_rates_ad_consistency(; ε::Float64=1e-4)
    r = Float64[0.1, 0.2, 0.3, 0.4, 0.5]
    fixed_rates(copy(r), tuple()) ≈ fixed_rates_ad(r, tuple()) || return false
    fe = ([1, 3, 4],)
    fixed_rates(copy(r), fe) ≈ fixed_rates_ad(r, fe) || return false
    fe2 = ([2, 4], [1, 3, 5])
    fixed_rates(copy(r), fe2) ≈ fixed_rates_ad(r, fe2) || return false

    transitions = ([1, 2], [2, 1])
    G = 2
    rtarget = [0.33, 0.19, 20.5, 1.0]
    nRNA = 60
    h = simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=20_000, nalleles=2)[1]
    data = RNAData{typeof(nRNA),typeof(h)}("", "", nRNA, h, 1.0, [])
    rinit = [0.1, 0.1, 0.1, 1.0]
    fittedparam = [1, 2, 3]
    model = load_model(
        data, rinit, prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0),
        fittedparam, tuple(), transitions, G, 0, 0, 0, "", 2, 10.0, Int[], rtarget[end], 0.02,
        prob_Gaussian, [], 1, tuple(), tuple(), nothing,
    )
    θ = Vector{Float64}(get_param(model))
    get_rates(θ, model) ≈ get_rates_ad(θ, model) || return false
    p1 = predictedfn(θ, data, model; steady_state_solver=:augmented, rates_fn=get_rates)
    p2 = predictedfn(θ, data, model; steady_state_solver=:augmented, rates_fn=get_rates_ad)
    p1 ≈ p2 || return false

    R = 1
    S = 1
    insertstep = 1
    rtarget2 = [0.02, 0.1, 0.5, 0.2, 0.1, 0.01]
    nhist = 24
    bins = collect(1:1.0:200.0)
    hs = simulator(
        rtarget2, transitions, G, R, S, insertstep;
        nalleles=2,
        nhist=nhist,
        totalsteps=80_000,
        bins=bins,
    )
    hRNA = div.(hs[1], 30)
    data2 = RNAOnOffData("test", "test", nhist, hRNA, bins, hs[2], hs[3], 1.0)
    rinit2 = fill(0.01, num_rates(transitions, R, S, insertstep))
    fittedparam2 = collect(1:length(rtarget2) - 1)
    model2 = load_model(
        data2,
        rinit2,
        prior_ratemean(transitions, R, S, insertstep, rtarget2[end], [], mean_elongationtime(rtarget2, transitions, R)),
        fittedparam2,
        tuple(),
        transitions,
        G,
        R,
        S,
        insertstep,
        "",
        2,
        10.0,
        Int[],
        rtarget2[end],
        0.05,
        prob_Gaussian,
        [],
        1,
        tuple(),
        tuple(),
        nothing,
    )
    θ2 = Vector{Float64}(get_param(model2))
    get_rates(θ2, model2) ≈ get_rates_ad(θ2, model2) || return false
    q1 = predictedfn(θ2, data2, model2; steady_state_solver=:augmented, rates_fn=get_rates)
    q2 = predictedfn(θ2, data2, model2; steady_state_solver=:augmented, rates_fn=get_rates_ad)
    q1 ≈ q2 || return false

    fsum(θv) = sum(get_rates_ad(θv, model))
    zg = Zygote.gradient(fsum, θ)[1]
    fdg = _test_central_grad(fsum, θ; ε=ε)
    isapprox(zg, fdg; rtol=1e-3, atol=1e-3) || return false

    return true
end

"""
    test_run_nuts_fit_smoke(; seed=42, n_samples=12, n_adapts=12)

Short finite-difference NUTS smoke test used by `test/runtests.jl`.
Returns `true` when the returned fit/statistics shapes match the expected `run_mh`-style layout.
"""
function test_run_nuts_fit_smoke(; seed::Int=42, n_samples::Int=12, n_adapts::Int=12)
    transitions = ([1, 2], [2, 1])
    G = 2
    rtarget = [0.33, 0.19, 20.5, 1.0]
    nRNA = 40
    h = simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=8_000, nalleles=2)[1]
    data = RNAData{typeof(nRNA),typeof(h)}("", "", nRNA, h, 1.0, [])
    rinit = [0.1, 0.1, 0.1, 1.0]
    fittedparam = [1, 2, 3]
    model = load_model(
        data, rinit, prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0),
        fittedparam, tuple(), transitions, G, 0, 0, 0, "", 2, 10.0, Int[], rtarget[end], 0.02,
        prob_Gaussian, [], 1, tuple(), tuple(), nothing,
    )
    opts = NUTSOptions(
        ;
        n_samples=n_samples,
        n_adapts=n_adapts,
        δ=0.8,
        gradient=:finite,
        fd_ε=1e-4,
        verbose=false,
        progress=false,
    )
    rng = MersenneTwister(seed)
    fits, stats, measures, nuts_info = run_nuts_fit(
        data, model, opts;
        rng=rng,
        steady_state_solver=:augmented,
    )
    size(fits.param, 2) == n_samples || return false
    length(fits.ll) == n_samples || return false
    fits.total == opts.n_adapts + n_samples || return false
    length(stats.meanparam) == size(fits.param, 1) || return false
    measures.waic isa Tuple || return false
    length(measures.ess) == size(fits.param, 1) || return false
    haskey(nuts_info, :nuts_stats) || return false
    haskey(nuts_info, :initial_θ) || return false

    out = run_NUTS(data, model, opts; rng=MersenneTwister(seed), steady_state_solver=:augmented)
    size(out[1].param, 2) == n_samples || return false
    return true
end

### end of functions used in runtest

# ════════════════════════════════════════════════════════════════════════════════════
# SECTION 4: EXTENDED & EXPERIMENTAL TESTS
# ════════════════════════════════════════════════════════════════════════════════════
# Variants, extensions, and experimental tests for specialized scenarios.
# Exported for interactive use but NOT called by Pkg.test().
# Useful for testing specific model configurations, hierarchical structures, and edge cases.

### developer-only benchmark helpers (not called by `runtests.jl`)

"""
    test_trace_subset_benchmark_keyword_bundle(; kwargs...)

Smoke-test the trace subset benchmark helpers with one shared keyword bundle,
including `fd_ε`. The finite-difference path uses `fd_ε`; the ForwardDiff and
Zygote subset benchmark helpers should accept it for API compatibility.

Use `scenario=:gr2r2` for a fast default smoke test or `scenario=:coupled_3x3`
to reproduce the uploaded coupled benchmark shape.

Returns a `NamedTuple` with `forwarddiff`, `finitediff`, and `zygote` results.
"""
function test_trace_subset_benchmark_keyword_bundle(
    ;
    scenario::Symbol=:gr2r2,
    seed::Int=2,
    totaltime::Float64=40.0,
    ntrials::Int=1,
    interval::Float64=5 / 3,
    fittedparam::AbstractVector{<:Integer}=Int[1, 2, 3],
    param_indices::AbstractVector{<:Integer}=Int[1],
    nruns::Int=1,
    warmup::Bool=false,
    steady_state_solver::Symbol=:augmented,
    hmm_checkpoint_steps::Union{Nothing,Integer}=8,
    fd_ε::Float64=1e-4,
)
    scen = if scenario === :gr2r2
        benchmark_inference_trace_gr2r2(
            ;
            seed=seed,
            totaltime=totaltime,
            ntrials=ntrials,
            interval=interval,
            fittedparam=collect(Int, fittedparam),
        )
    elseif scenario === :coupled_3x3
        benchmark_inference_trace_coupled_3x3(
            ;
            seed=seed,
            totaltime=totaltime,
            ntrials=ntrials,
            interval=interval,
            fittedparam=collect(Int, fittedparam),
        )
    else
        throw(ArgumentError("scenario must be :gr2r2 or :coupled_3x3 (got $(scenario))"))
    end
    kw = (
        ;
        param_indices=collect(Int, param_indices),
        nruns=nruns,
        warmup=warmup,
        steady_state_solver=steady_state_solver,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
        fd_ε=fd_ε,
    )
    forwarddiff = benchmark_trace_forwarddiff_gradient(scen; kw...)
    finitediff = benchmark_trace_finitediff_gradient(scen; kw...)
    zygote = benchmark_trace_zygote_subset_gradient(scen; kw...)
    return (; forwarddiff, finitediff, zygote)
end
=#

### functions to be used in the future
"""
    test_fit_tracejoint_hierarchical(; coupling, G, R, S, insertstep, transitions, rtarget, rinit, hierarchical, method, nsamples, onstates, totaltime, ntrials, fittedparam, propcv, cv, interval, noisepriors, maxtime, decayrate)

Fit a simulated joint trace dataset for coupled models using a hierarchical model and compare to the target.

# Arguments
- `coupling`, `G`, `R`, `S`, `insertstep`, `transitions`: Model structure and coupling.
- `rtarget`, `rinit`: Rate parameters.
- `hierarchical`, `method`: Hierarchical model and method options.
- `nsamples`, `onstates`, `totaltime`, `ntrials`: Data and simulation options.
- `fittedparam`, `propcv`, `cv`, `interval`, `noisepriors`, `maxtime`, `decayrate`: Fitting options.

# Returns
- Tuple of (fitted rates, target rates, fits, stats, measures, model, data, options).
"""
function test_fit_tracejoint_hierarchical(; coupling=((1, 2), [(1, 2, 2, 1)], [:free]), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])),   rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 0., 0.05, 1., 0.05, 0.03, 0.1, 0.5, 0.2, 0.1, 0., 0.05, 1., 0.05, -0.5], rinit=Float64[], hierarchical=(2, [8, 17], tuple()), method=(Tsit5(), true), nsamples=20000, onstates=Int[], totaltime=1000.0, ntrials=20, fittedparam=Int[1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 21], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([0., .1, 1., 0.05], [0., 0.1, 1., 0.1]), maxtime=300.0, decayrate=1.0)
    rh = 50.0 .+ 10 * randn(ntrials)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials; hierarchical=(6, rh), noiseparams=[4, 4])
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], [0.0, 0.0], 1), [1, 2])
    priormean = set_priormean([], transitions, R, S, insertstep, decayrate, noisepriors, mean_elongationtime(rtarget, transitions, R), hierarchical, coupling, nothing)
    rinit = isempty(hierarchical) ? set_rinit(rinit, priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), coupling, nothing; nhypersets=hierarchical[1])
    fittedparam = set_fittedparam(fittedparam, "tracejoint", transitions, R, S, insertstep, noisepriors, coupling, nothing)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, hierarchical, coupling, nothing)
    options = StochasticGene.MHOptions(nsamples, 0, maxtime_seconds(maxtime), 1.0)
    fits, stats, measures = run_mh(data, model, options)
    rfit = StochasticGene.get_rates(fits.parml, model)
    return rfit[1:length(rtarget)], rtarget, fits, stats, measures, model, data, options
end

"""
    test_fit_tracejoint_3unit(; coupling, G, R, S, insertstep, transitions, rtarget, rinit, nsamples, onstates, totaltime, ntrials, fittedparam, propcv, cv, interval, noisepriors, maxtime, method)

Fit a simulated joint trace dataset for 3 coupled units: units 1 and 2 same as test_fit_tracejoint (two G=2, R=2/1 models);
unit 3 is a G=3 telegraph (R=0) that inhibits units 1 and 2 and has no trace (hidden). Only units 1 and 2 have observations;
`units=[1, 2]` is passed to the simulator and stored in data.units.

# Arguments
- `coupling`: (unit_model, connections, sign_modes) with unit_model=(1,2,3), connections=[(3,1,1,2), (3,3,2,2)], sign_modes=[:inhibit,:inhibit].
- `G=(2, 2, 3)`, `R=(2, 1, 0)`, `S=(1, 0, 0)`, `insertstep=(1, 1, 0)`.
- `transitions`: first two = ([1,2],[2,1]); third = ([1,2],[2,1],[2,3],[3,2],[1,3],[3,1]).
- `units`: unit indices that have traces (default [1, 2]; unit 3 is hidden).
- `fixedeffects`: tuple of tied indices: unit 3 rates (all tied to first) and the two coupling constants (tied). Default builds ([unit3_rate_start, ..., unit3_rate_end], [coupling1, coupling2]) so all model-3 rates share one parameter (0.1) and both coupling params are equal.
- `fittedparam`: indices to fit; default is `[unit3_rate_start, coupling_start]` (unit-3 master rate and coupling master).
- `rtarget`, `rinit`, `nsamples`, `totaltime`, `ntrials`, `propcv`, `interval`, `noisepriors`, `maxtime`, `method`: as in test_fit_tracejoint.

# Returns
- Tuple of (lower, rtarget[fittedparam], upper) or (fits, stats, measures, model, data, options) for inspection.
"""
function test_fit_tracejoint_3unit(;
    coupling=((1, 2, 3), [(3, 1, 1, 2), (3, 3, 2, 2)], [:inhibit, :inhibit]),
    G=(2, 2, 3),
    R=(2, 1, 0),
    S=(1, 0, 0),
    insertstep=(1, 1, 1),
    transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1]), ([1, 2], [2, 1], [2, 3], [3, 2], [1, 3], [3, 1])),
    units=[1, 2],
    fixedeffects=nothing,
    rtarget=nothing,
    rinit=Float64[],
    nsamples=5000,
    onstates=Int[],
    totaltime=1000.0,
    ntrials=50,
    fittedparam=Int[],
    propcv=0.05,
    cv=10.0,
    interval=1.0,
    noisepriors=([0.0, 0.1, 1.0, 0.1], [0.0, 0.1, 1.0, 0.1], [0.0, 0.1, 1.0, 0.1]),
    maxtime="60m",
    method=Tsit5(),
)
    nrates_per_unit = StochasticGene.num_rates(transitions, R, S, insertstep)
    ncoupling = StochasticGene.ncoupling(coupling)
    n_noise = sum(length.(noisepriors))
    n_total = sum(nrates_per_unit) + n_noise + ncoupling
    # Full vector layout: unit1 rates (1..n1), unit1 noise, unit2 rates, unit2 noise, unit3 rates, unit3 noise, coupling
    unit3_rate_start = sum(nrates_per_unit[1:2]) + sum(length.(noisepriors)[1:2]) + 1
    unit3_rate_end = unit3_rate_start + nrates_per_unit[3] - 1
    coupling_start = sum(nrates_per_unit) + n_noise + 1
    if fixedeffects === nothing
        fixedeffects = (collect(unit3_rate_start:unit3_rate_end), collect(coupling_start:coupling_start + ncoupling - 1))
    end
    if rtarget === nothing
        rtarget = fill(0.1, n_total)
        rtarget[sum(nrates_per_unit[1:1])] = 1.0
        rtarget[sum(nrates_per_unit[1:2])] = 1.0
        rtarget[sum(nrates_per_unit[1:3])] = 1.0
        rtarget[unit3_rate_start:unit3_rate_end] .= 0.1
        for k in 1:ncoupling
            rtarget[coupling_start - 1 + k] = -0.2
        end
    end
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials; observed_units=units, noiseparams=[4, 4, 0])
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], fill(0.0, length(units)), 1), units)
    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors, mean_elongationtime(rtarget, transitions, R), tuple(), coupling, nothing)
    rinit = isempty(rinit) ? copy(rtarget) : rinit
    if isempty(fittedparam)
        # Fit unit-3 (master) rate and (master) coupling parameter; other unit-3 rates and second coupling are tied via fixedeffects
        fittedparam = [unit3_rate_start, coupling_start]
    end
    # Decay index per unit in block layout: unit i's rates are at block_starts[i] : block_starts[i]+nrates_per_unit[i]-1
    block_ends = cumsum([nrates_per_unit[i] + length(noisepriors[i]) for i in eachindex(R)])
    block_starts = vcat(1, block_ends[1:end-1] .+ 1)
    decayrate = tuple((rtarget[block_starts[i] + nrates_per_unit[i] - 1] for i in eachindex(R))...)
    model = load_model(data, rinit, priormean, fittedparam, fixedeffects, transitions, G, R, S, insertstep, "", 1, 10.0, Int[], decayrate, propcv, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing)
    options = StochasticGene.MHOptions(nsamples, 100, maxtime_seconds(maxtime), 1.0)
    fits, stats, measures = run_mh(data, model, options)
    println(fits.accept," ",fits.total)
    lower, upper = compute_quantile_bounds(stats.qparam[1, :], stats.qparam[3, :], fits.accept, fits.total)
    return lower, rtarget[fittedparam], upper
end

"""
    test_fit_trace_grid(; grid, G, R, S, insertstep, transitions, rtarget, rinit, nsamples, onstates, totaltime, ntrials, fittedparam, propcv, cv, interval, weight, nframes, noisepriors, maxtime)

Fit a simulated trace grid dataset using the provided parameters and compare to the target.

# Arguments
- `grid`, `G`, `R`, `S`, `insertstep`, `transitions`: Model structure.
- `rtarget`, `rinit`: Rate parameters.
- `nsamples`, `onstates`, `totaltime`, `ntrials`: Data and simulation options.
- `fittedparam`, `propcv`, `cv`, `interval`, `weight`, `nframes`, `noisepriors`, `maxtime`: Fitting options.

# Returns
- Tuple of (fits, stats, measures, data, model, options).
"""
function test_fit_trace_grid(; grid=4, G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 20, 0.2], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [50, 15, 200, 20]; 0.1], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[1, 2, 3, 4, 5, 6, 7, 8, 11], propcv=0.01, cv=100.0, interval=1.0, weight=0, nframes=1, noisepriors=[50, 15, 200, 70], maxtime=10.0)
    traces = sim_grid(r=rtarget[1:end-1], p=rtarget[end], Ngrid=grid, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, totaltime=totaltime, interval=interval, ntrials=ntrials)
    data = StochasticGene.TraceData("tracegrid", "test", interval, (traces, [], weight, nframes), Int[])
    priormean = set_priormean([], transitions, R, S, insertstep, rtarget[num_rates(transitions, R, S, insertstep)], noisepriors, mean_elongationtime(rtarget, transitions, R), tuple(), tuple(), grid)
    rinit = set_rinit(rinit, priormean)
    fittedparam = set_fittedparam(fittedparam, "tracegrid", transitions, R, S, insertstep, noisepriors, tuple(), grid)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian_grid, noisepriors, Tsit5(), tuple(), tuple(), grid)
    options = StochasticGene.MHOptions(nsamples, 0, maxtime_seconds(maxtime), 1.0)
    fits, stats, measures = run_mh(data, model, options)
    fits, stats, measures, data, model, options
end

"""
    test_fit_trace_grid_hierarchical(; grid, G, R, S, insertstep, transitions, rtarget, rinit, nsamples, onstates, hierarchical, totaltime, ntrials, fittedparam, propcv, cv, interval, weight, nframes, noisepriors, maxtime)

Fit a simulated trace grid dataset using a hierarchical model and compare to the target.

# Arguments
- `grid`, `G`, `R`, `S`, `insertstep`, `transitions`: Model structure.
- `rtarget`, `rinit`: Rate parameters.
- `hierarchical`: Hierarchical model options.
- `nsamples`, `onstates`, `totaltime`, `ntrials`: Data and simulation options.
- `fittedparam`, `propcv`, `cv`, `interval`, `weight`, `nframes`, `noisepriors`, `maxtime`: Fitting options.

# Returns
- Tuple of (fits, stats, measures, data, model, options).
"""
function test_fit_trace_grid_hierarchical(; grid=4, G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 20, 0.2], rinit=[], nsamples=5000, onstates=Int[], hierarchical=(2, [7], tuple()), totaltime=1000.0, ntrials=10, fittedparam=[1, 2, 3, 4, 5, 11], propcv=0.01, cv=100.0, interval=1.0, weight=0, nframes=1, noisepriors=[50, 15, 200, 70], maxtime=10.0)
    traces = sim_grid(r=rtarget[1:end-1], p=rtarget[end], Ngrid=grid, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, totaltime=totaltime, interval=interval, ntrials=ntrials)
    data = StochasticGene.TraceData("tracegrid", "test", interval, (traces, [], weight, nframes), Int[])
    priormean = set_priormean([], transitions, R, S, insertstep, rtarget[num_rates(transitions, R, S, insertstep)], noisepriors, mean_elongationtime(rtarget, transitions, R), hierarchical, tuple(), grid)
    rinit = isempty(hierarchical) ? set_rinit(rinit, priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), tuple(), grid; nhypersets=hierarchical[1])
    fittedparam = set_fittedparam(fittedparam, "tracegrid", transitions, R, S, insertstep, noisepriors, tuple(), grid)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian_grid, noisepriors, (Tsit5(), true), hierarchical, tuple(), grid)
    options = StochasticGene.MHOptions(nsamples, 0, maxtime_seconds(maxtime), 1.0)
    fits, stats, measures = run_mh(data, model, options)
    fits, stats, measures, data, model, options
end

### under development

function test_fit_trace_forced(; datapath="data/forced/2", label="trace-HBEC-nstate", gene="MYC", datacond=["enhancer", "gene"], coupling=((1, 2), [(1, 2, 2, 1)]), traceinfo=(1.0, 1.0, -1, 1.0, 0.5), G=2, R=2, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.05, 0.2, 0.1, 0.15, 0.1, 1.0, 50, 5, 50, 5], nsamples=10000000, onstates=Int[], totaltime=1000.0, ntrials=100, fittedparam=[1:num_rates(transitions, R, S, insertstep)-1; num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+1], propcv=0.01, cv=100.0, noisepriors=[0.0, 0.1, 1.0, 0.1], nchains=1, zeromedian=[false, true], maxtime=10.0, initprior=0.1)

    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors, mean_elongationtime(rtarget, transitions, R), tuple(), tuple(), nothing)
    elongationtime = mean_elongationtime(rtarget, transitions, R)
    priormean = [fill(0.01, length(transitions)); initprior; fill(R / elongationtime, R); fill(0.05, max(0, S - insertstep + 1)); 1.0; noisepriors; 0.01]
    priorcv = [fill(1.0, length(transitions)); 0.1; fill(0.1, R); fill(0.1, max(0, S - insertstep + 1)); 1.0; [0.5, 0.5, 0.1, 0.1]; 1.]
    rinit = isempty(tuple()) ? set_rinit([], priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, maxtime_seconds(maxtime), 1.0)

    data = load_data_trace(datapath, label, gene, datacond, traceinfo, :tracejoint, 1, zeromedian)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, priorcv, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, Tsit5(), tuple(), coupling, nothing, zeromedian)
    data2 = load_data_trace(datapath, label, gene, datacond[2], traceinfo, :trace, 1, zeromedian[2])
    model2 = load_model(data2, rinit[1:10], priormean[1:10], fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, priorcv[1:10], Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, Tsit5(), tuple(), tuple(), nothing, zeromedian)

    ll2 = loglikelihood(get_param(model2), data2, model2)
    ll1 = loglikelihood(get_param(model), data, model)
    # return ll1, ll2
    fits, stats, measures = run_mh(data, model, options, nchains)
    fits2, stats2, measures2 = run_mh(data2, model2, options, nchains)
    # nrates = num_rates(model)
    # stats.medparam[1:nrates-1], rtarget[1:nrates-1]
    return fits, stats, measures, data, model, options, fits2, stats2, measures2, data2, model2
end

function sim_tracejoint(; coupling=((1, 2), [(1, 2, 2, 1)]), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 1.0, 0., 0.05, 1.0, 0.05, 0.05, 0.023, 0.15, 0.2, 1.0, 0.0, 0.05, 1.0, 0.05, 0.], rinit=Float64[], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=Int[21], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([0., 0.1, 1., 0.1], [0., 0.1, 1., 0.1]), maxtime="60m", method=Tsit5())
    simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials)
end

function test_tracejoint(traces; coupling=((1, 2), [(1, 2, 2, 1)]), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 1.0, 0., 0.05, 1.0, 0.05, 0.03, 0.1, 0.5, 0.2, 1.0, 0.0, 0.05, 1.0, 0.05, 0.], fittedparam=Int[21], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([0., 0.1, 1., 0.1], [0., 0.1, 1., 0.1]), method=Tsit5())

    data = StochasticGene.TraceData("tracejoint", "test", interval, (traces, [], [0.0, 0.0], 1), Int[])
    fittedparam = set_fittedparam(fittedparam, "trace", transitions, R, S, insertstep, noisepriors, coupling, nothing)
    model = load_model(data, rtarget, rtarget, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], 1.0, propcv, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing, [true, true])
    ll = loglikelihood(get_param(model), data, model)

    n = num_rates(transitions, R, S, insertstep)

    trace1 = [t[:,1] for t in traces]
    data1 = StochasticGene.TraceData("trace", "test", interval, (trace1, [], 0., 1), Int[])
    r1 = rtarget[1:n[1]+4]
    model1 = load_model(data1, r1, r1, [1,2], tuple(), transitions[1], G[1], R[1], S[1], insertstep[1], "", 1, 10.0, Int[], 1., propcv, prob_Gaussian, noisepriors[1], Tsit5(), tuple(), tuple(), nothing, true)
    ll1 = loglikelihood(get_param(model1), data1, model1)

    trace2 = [t[:, 2] for t in traces]
    data2 = StochasticGene.TraceData("trace", "test", interval, (trace2, [], 0., 1), Int[])
    r2 = rtarget[n[1]+5:n[1]+n[2]+9]
    model2 = load_model(data2, r2, r2, [1,2], tuple(), transitions[2], G[2], R[2], S[2], insertstep[2], "", 1, 10.0, Int[], 1., propcv, prob_Gaussian, noisepriors[2], Tsit5(), tuple(), tuple(), nothing, true)
    ll2 = loglikelihood(get_param(model2), data2, model2)

    return ll[1], ll1[1]+ll2[1], ll1[1], ll2[1]

end

function test_tracejoint_hierarchical(; coupling=((1, 2), [(1, 2, 2, 1)]), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 0., 0.05, 1., 0.05, 0.03, 0.1, 0.5, 0.2, 0.1, 0., 0.05, 1., 0.05, 0.], rinit=Float64[], hierarchical=(2, [8, 17], tuple()), method=(Tsit5(), true), nsamples=20000, onstates=Int[], totaltime=1000.0, ntrials=20, fittedparam=Int[1, 2, 3, 4, 5, 6], interval=1.0, noisepriors=([0., 0.1, 1., 0.05], [0., 0.1, 1., 0.1]), decayrate=1.0)
    # trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials)
    rh = 50.0 .+ 10 * randn(ntrials)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials, hierarchical=(6, rh))

    trace1 = [t[:,1] for t in trace]
    trace2 = [t[:,2] for t in trace]

    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [0.0, 0.0], 1), Int[])
    data1 = StochasticGene.TraceData("trace","test",interval, (trace1,[],0., 1), Int[])
    data2 = StochasticGene.TraceData("trace","test",interval, (trace2,[],0., 1), Int[])

    hierarchical = tuple()
    model = load_model(data, rtarget, rtarget, [1], tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], 1., .01, prob_Gaussian, noisepriors, method, hierarchical, coupling, nothing)

    ll = loglikelihood(get_param(model), data, model)

    return ll
    n = num_rates(transitions, R, S, insertstep)
    r1 = rinit[1:n[1]+4]
    r2 = rinit[n[1] + 5: n[1] + n[2] + 9]
    priormean1 = priormean[1:n[1]+4]
    priormean2 = priormean[n[1] + 5: n[1] + n[2] + 9]
    print(r1)
    print(r2)

    data1 = StochasticGene.TraceData("trace","test",interval, (trace,[],0., 1), Int[])
    data2 = StochasticGene.TraceData("trace","test",interval, (trace,[],0., 1), Int[])
    model1 = load_model(data1, r1, priormean1, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rinit[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, hierarchical, coupling, nothing)
    model2 = load_model(data2, r2, priormean2, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rinit[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, hierarchical, coupling, nothing)
    ll1 = loglikelihood(get_param(model1), data1, model1)
    ll2 = loglikelihood(get_param(model2), data2, model2)
    return ll, ll1, ll2
end

### development test functions

test_fit0(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=true) = fit(nchains, "trace", String[], "data/inhibition/control/", "MYC", "HBEC", "gene", (1.6666666666666667, 1.0, -1, 0.92, 0.0), "test", "test", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3, 4, 5, 6, 8, 9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 0.0, 0.2, 1.0, 0.2], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.01, 0.01], 1, Int64[], 1.0, "", prob_Gaussian, [0.0, 0.2, 1.0, 0.2], (), "ml", propcv, 200000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5(), zeromedian)

test_fit1(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=false, warmup=0) = fit(nchains, "trace", String[], "data/inhibition/control/", "MYC", "HBEC", "gene", (1.6666666666666667, 1.0, -1, 0.92, 0.5), "test1", "test1", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3, 4, 5, 6, 8, 9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 0.5, 0.25, 1.0, 0.2], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.5, 0.5, 0.1, 0.1], 1, Int64[], 1.0, "", prob_Gaussian, [0.5, 0.2, 1.0, 0.2], (), "ml", propcv, 20000, warmup, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5(), zeromedian)

test_fit2(; maxtime=6.0, propcv=0.02, zeromedian=true) = fit(1, "trace", String[], "data/inhibition/inhibition/", "MYC", "HBEC", "gene", (1.6666666666666667, 1.0, -1, 0.78, 0.0), "inhibition-test", "inhibition-test", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3, 4, 5, 6, 7, 9], (), ([1, 2], [2, 1]), 2, 4, 0, 3, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.8, 0.8, 0.8, 0.8, 1.0, 0.0, 0.25, 1.0, 0.25], [1.0, 1.0, 0.2, 0.1, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 0.1], 1, Int64[], 1.0, "", prob_Gaussian, [0.0, 0.25, 1.0, 0.25], (), "ml", propcv, 2000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5(), zeromedian, 3, 1);

test_fitt(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=false) = fit(nchains, "trace", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", "gene", (1.0, 1.0, -1, 1.0, 0.5), "genetogther", "genetogther", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3, 4, 5, 6, 8, 9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 0.5, 0.2, 1.0, 0.2], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.5, 0.5, 0.1, 0.1], 1, Int64[], 1.0, "", prob_Gaussian, [0.5, 0.2, 1.0, 0.2], (), "ml", propcv, 200000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5(), zeromedian)

test_fit_h(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=true, warmup=0) = fit(nchains, "trace", String[], "data/inhibition/control/", "MYC", "HBEC", "gene", (1.0, 1.0, -1, 0.92, 0.0), "testh", "testh", "trace-HBEC-nstate-h_gene", "trace-HBEC-nstate-h_gene", [1, 2, 3, 4, 5, 6], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 0.0, 0.2, 1.0, 0.1, 1.0, 1.0, 1.0, 0.1, 0.1, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 0.1, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0], 1, Int64[], 1.0, "", prob_Gaussian, [0.0, 0.2, 1.0, 0.1], (2, [8], ()), "ml", propcv, 100000, warmup, 0, 1.0, 100.0, 1.0, false, false, false, (Tsit5(), true), zeromedian)

test_fit_h2(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=true) = fit(nchains, "trace", String[], "data/inhibition/inhibition/", "MYC", "HBEC", "gene", (1.0, 1.0, -1, 0.78, 0.0), "inhibition-test", "inhibition-test", "trace-HBEC-nstate-h_gene", "trace-HBEC-nstate-h_gene", [1, 2, 3, 4, 5, 6, 8, 9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 0.0, 0.25, 1.0, 0.25, 1.0, 1.0, 1.0, 0.0, 0.1, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 0.1, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0], 1, Int64[], 1.0, "", prob_Gaussian, [0.0, 0.25, 1.0, 0.25], (2, [8], ()), "ml", propcv, 200000, 0, 0, 1.0, 100.0, 1.0, false, false, false, (Tsit5(), true), zeromedian)

test_fit_coupleda(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=true) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [0.0, 0.0]), "test", "test", "tracejoint-HBEC-nstate_enhancer-geneR5", "tracejoint-HBEC-nstate_enhancer-geneR5", [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 28], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (3, 3), (1, 0), (1, 1), ((1, 2), [(1, 4, 2, 5), (1, 5, 2, 5), (1, 6, 2, 5)]), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([0.0, 0.2, 1.0, 0.1], [0.0, 0.2, 1.0, 0.1]), (), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5(), zeromedian)

test_fit_coupled2a(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=false) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [0.5, 0.5]), "test1", "test1", "tracejoint-HBEC-nstate-h_enhancer-gene11", "tracejoint-HBEC-nstate-h_enhancer-gene11", [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 28], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (3, 3), (1, 0), (1, 1), ((1, 2), [(1, 1, 2, 1)]), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([0.5, 0.2, 1.0, 0.1], [0.5, 0.2, 1.0, 0.1]), (), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, true, false, false, Tsit5(), zeromedian)

test_fit_coupled(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=false) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [0.5, 0.5]), "test1", "test1", "tracejoint-HBEC-nstate_enhancer-gene11", "tracejoint-HBEC-nstate-h_enhancer-gene11", [1, 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18, 24], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (1, 1), (1, 0), (1, 1), ((1, 2), [(1, 1, 2, 1)]), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([0.5, 0.2, 1.0, 0.1], [0.5, 0.2, 1.0, 0.1]), (), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5(), zeromedian)
test_fit_coupled0(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=true) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [0.0, 0.0]), "test", "test", "tracejoint-HBEC-nstate_enhancer-gene11", "tracejoint-HBEC-nstate-h_enhancer-gene11", [1, 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18, 24], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (1, 1), (1, 0), (1, 1), ((1, 2), [(1, 1, 2, 1)]), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([0.0, 0.2, 1.0, 0.1], [0.0, 0.2, 1.0, 0.1]), (), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5(), zeromedian)

test_fit_coupled_h(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=false) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [0.5, 0.5]), "test1", "test1", "tracejoint-HBEC-nstate-h_enhancer-gene11", "tracejoint-HBEC-nstate-h_enhancer-gene11", [1, 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18, 24], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (1, 1), (1, 0), (1, 1), ((1, 2), [(1, 1, 2, 1)]), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([0.5, 0.2, 1.0, 0.1], [0.5, 0.2, 1.0, 0.1]), (2, [9, 20], ()), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, (Tsit5(), true), zeromedian)
test_fit_coupled_h0(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=true) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [0.0, 0.0]), "test", "test", "tracejoint-HBEC-nstate-h_enhancer-gene11", "tracejoint-HBEC-nstate-h_enhancer-gene11", [1, 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18, 24], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (1, 1), (1, 0), (1, 1), ((1, 2), [(1, 1, 2, 1)]), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([0.0, 0.2, 1.0, 0.1], [0.0, 0.2, 1.0, 0.1]), (2, [9, 20], ()), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, (Tsit5(), true), zeromedian)

test_rnacount(; nchains=1, gene="RPLP1", datapath="data/U3AS4/counts/WT-UTR", maxtime="60m", propcv=0.01) = fit(nchains=nchains, datatype="rnacount", transitions=([1, 2], [2, 1]), G=2, R=0, S=0, insertstep=1, datapath=datapath, cell="U3A", gene=gene, datacond="WT-UTR", infolder="test", resultfolder="test", fittedparam=[1, 2, 3], maxtime=maxtime, ratetype="ml", propcv=propcv);

test_rna(; nchains=1, gene="RPLP1", datapath="data/U3AS4/histograms/", maxtime="60m", propcv=0.01) = fit(nchains=nchains, datatype="rna", transitions=([1, 2], [2, 1]), G=2, R=0, S=0, insertstep=1, datapath=datapath, cell="U3A", gene=gene, datacond="WT-UTR", infolder="test", resultfolder="test", fittedparam=[1, 2, 3], maxtime=maxtime, ratetype="ml", propcv=propcv);

test_fit_G(; nchains=1, maxtime="60m", propcv=0.01) = fit(nchains, "trace", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", "gene", (1.0, 1.0, -1, 0.86, [25000.0, 0.0]), "genetogther", "genetogther", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3], (), ([1, 2], [2, 1]), 2, 0, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 1.0, 20000.0, 15000.0, 80000.0, 4000.0, 1.0, 1.0, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25], [5.0, 5.0, 0.2, 1.0, 0.5, 0.5, 1.0, 0.1, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0], 1, [2], 1.0, "", prob_Gaussian, [20000.0, 15000.0, 80000.0, 4000.0], (), "ml", propcv, 200000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5())

fit_myc3(; nchains=1, maxtime="60m", propcv=0.05) = fit(nchains=nchains, datatype="rna", transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=0, S=0, insertstep=1, datapath="data/FISH", cell="HBEC", gene="MYC", datacond="gene", infolder="test", resultfolder="test", fittedparam=[1, 2, 3, 4, 5], maxtime=maxtime, ratetype="ml", propcv=propcv, burst=true);

fit_myc2(; nchains=1, maxtime="60m", propcv=0.05) = fit(nchains=nchains, datatype="rna", transitions=([1, 2], [2, 1]), G=2, R=0, S=0, insertstep=1, datapath="data/FISH", cell="HBEC", gene="MYC", datacond="gene", infolder="test", resultfolder="test", fittedparam=[1, 2, 3], maxtime=maxtime, ratetype="ml", propcv=propcv, burst=true);

# @time fit(16, "trace", String[], "data/inhibition/control/", "MYC", "HBEC", "gene", (1.6666666666666667, 1.0, -1, 0.92, [290.0, 140.0]), "control-2025-03-15", "control-2025-03-15", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene",[1, 2,3, 4, 5, 6, 8, 9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", 43000.0, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 290.0, 140.0, 1200.0, 175.0], [2.0, 2.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.5, 0.5, 0.1, 0.1], 1, Int64[], 1.0, "", prob_Gaussian, [290.0, 140.0, 1200.0, 175.0], (), "ml", 0.02, 2000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5())


##### Experimental functions

function aic_onstates(ratefile, datapath, gene, datacond, traceinfo, label, fittedparam, transitions, G, R, S, insertstep, hierarchical, ratetype)
    r = readrates(ratefile, ratetype)
    data = load_data_trace(datapath, label, gene, datacond, traceinfo, dt, datacol, zeromedian)
    model = load_model(data, r, r, fittedparam, [], transitions, G, R, S, insertstep, "", 1, 0.1, Int[], 1., 0.1, prob_Gaussian, [0.0, 0.2, 1.0, 0.2], Tsit5(), hierarchical, tuple(), nothing, zeromedian, 1)
    aic_onstates(r, data, model)
end

function aic_onstates(param, data, model::AbstractGRSMmodel)
    r = prepare_rates(param, model)
    components = get_components(model, data)
    aic_onstates(r, components, model.reporter, data.interval, data.trace, model.method)
end

function aic_onstates(r::Tuple{T1,T2}, components::TComponents, reporter::HMMReporter, interval, trace, method=Tsit5()) where {T1,T2}
    rates, noiseparams = r
    a, p0 = make_ap(rates, interval, components, method)
    onstates = reporter.per_state .> 0
    ll_on, ll = _ll_onstates(a, p0, set_d(noiseparams, reporter), trace[1], onstates)
    2 * length(rates) - 2 * ll_on, 2 * length(rates) - 2 * ll, ll_on, ll
end

function aic_onstates(r::Tuple{T1,T2,T3,T4,T5,T6}, components::TComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method::Tuple=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6}
    rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper = r
    a, p0 = make_ap(rshared[1], interval, components, method[1])
    onstates = reporter.per_state .> 0
    if method[2]
        ll_on, ll = _ll_onstates(noiseindividual, a, p0, reporter, trace[1], onstates)
        aic = 2 * length(rshared[1]) - 2 * ll_on
    else
        # ll_on = _ll_on(rindividual, noiseindividual, interval, components, reporter, trace[1], method[1], onstates)
        # aic = 2 * length(rshared) * (length(rindividual) + 1) - 2 * ll_on
    end
    return aic, ll
end

function _ll_onstates(a::Matrix, p0::Vector, d, traces, onstates)
    ll_on = 0
    ll = 0
    for i in eachindex(traces)
        b = set_b(traces[i], d)
        α, C = forward(a, b, p0)
        # Normalize by dividing by C to get proper likelihood
        ll_on += log(sum(max.(α[onstates, end], 0.))) - sum(log.(C))
        ll -= sum(log.(C))
    end
    ll_on, ll
end

function _ll_onstates(noiseparams::Vector, a::Matrix, p0::Vector, reporter, traces, onstates)
    ll_on = 0
    ll = 0
    for i in eachindex(traces)
        d = set_d(noiseparams[i], reporter)
        b = set_b(traces[i], d)
        α, C = forward(a, b, p0)
        ll_on += log(sum(max.(α[onstates, end], 0.)))
        ll -= sum(log.(C))
    end
    ll_on, ll
end


###### test autodiff


function test_trace(; traceinfo=(1.0, 1.0, -1, 1.0, 0.5), G=2, R=2, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.05, 0.2, 0.1, 0.15, 0.1, 1.0, 50, 5, 50, 5], nsamples=10000000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[1:num_rates(transitions, R, S, insertstep)-1; num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+1], propcv=0.01, cv=100.0, noisepriors=[0.0, 0.1, 1.0, 0.1], zeromedian=true, maxtime=100.0, initprior=0.1)
    tracer = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, traceinfo[1], totaltime, ntrials)
    trace, tracescale = zero_median(tracer, zeromedian)
    nframes = round(Int, mean(size.(trace, 1)))  #mean number of frames of all traces
    if length(traceinfo) > 3 && traceinfo[4] != 1.0
        weight = set_trace_weight(traceinfo)
        background = set_trace_background(traceinfo)
    else
        weight = 0.0
        background = 0.0
    end
    data = TraceData{String,String,Tuple}("trace", "gene", traceinfo[1], (trace, background, weight, nframes, tracescale), Int[])
    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors, mean_elongationtime(rtarget, transitions, R), tuple(), tuple(), nothing)
    elongationtime = mean_elongationtime(rtarget, transitions, R)
    priormean = [fill(0.01, length(transitions)); initprior; fill(R / elongationtime, R); fill(0.05, max(0, S - insertstep + 1)); 1.0; noisepriors]
    priorcv = [fill(1.0, length(transitions)); 0.1; fill(0.1, R); fill(0.1, max(0, S - insertstep + 1)); 1.0; [0.5, 0.5, 0.1, 0.1]]

    rinit = isempty(tuple()) ? set_rinit([], priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), tuple(), nothing)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, priorcv, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, Tsit5(), tuple(), tuple(), nothing, zeromedian)
    return data, model
end

# function test_ad_loglikelihood()
#     # Use your existing test function to get data/model/param
#     data, model = test_trace()  # or whatever function returns these
#     param = get_param(model)
#     ll_wrap = param -> loglikelihood(param, data, model)
#     grad = Zygote.gradient(ll_wrap, param)
#     return grad
# end

# function ll_wrap(param)
#     loglikelihood(param, data, model)
# end

# grads = Zygote.gradient(ll_wrap, param)

##### Experimental / exploratory helpers (scrapyard; not used by `runtests.jl`) #####




function test_tracejoint_h(; coupling=((1, 2), [(1, 2, 2, 1)]), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 30, 100, 20, 0.03, 0.1, 0.5, 0.2, 0.1, 50, 30, 100, 20, -0.5], rinit=Float64[], hierarchical=(2, [8, 17], tuple()), method=(Tsit5(), true), nsamples=20000, onstates=Int[], totaltime=1000.0, ntrials=20, fittedparam=Int[1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 21], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([50, 30, 100, 20], [50, 30, 100, 20]), maxtime=300.0, decayrate=1.0)
    # trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials)
    rh = 50.0 .+ 10 * randn(ntrials)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials, hierarchical=(6, rh))
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], [0.0, 0.0], 1), Int[])
    priormean = set_priormean([], transitions, R, S, insertstep, decayrate, noisepriors, mean_elongationtime(rtarget, transitions, R), hierarchical, coupling, nothing)
    rinit = isempty(hierarchical) ? set_rinit(rinit, priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), coupling, nothing; nhypersets=hierarchical[1])
    fittedparam = set_fittedparam(fittedparam, "tracejoint", transitions, R, S, insertstep, noisepriors, coupling, nothing)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, hierarchical, coupling, nothing)
    options = StochasticGene.MHOptions(nsamples, 0, maxtime_seconds(maxtime), 1.0)
    return data, model, options
end


function test_tracejoint(; coupling=((1, 2), [(1, 2, 2, 1)]), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 30, 100, 20, 0.02, 0.05, 0.2, 0.2, 0.1, 50, 30, 100, 20, -0.5], rinit=Float64[], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=Int[], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([100, 50, 200, 100], [100, 50, 200, 100]), maxtime=300.0, method=Tsit5())
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials; verbose=false)
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], 0.0, 1), Int[])
    rm = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, [5.0, 5.0], coupling)
    rinit = set_rinit(rinit, rm)
    fittedparam = StochasticGene.set_fittedparam(fittedparam, data.label, transitions, R, S, insertstep, noisepriors, coupling, nothing)
    model = load_model(data, rinit, rm, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing)
    options = StochasticGene.MHOptions(nsamples, 0, maxtime_seconds(maxtime), 1.0)
    return data, model, options
end

function test_predicted_states(; coupling=((1, 2), [(1, 2, 2, 1)]), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.1, 0.3, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 20, 200, 20, 0.1, 0.2, 0.4, 0.2, 0.1, 50, 20, 200, 20, -0.5], rinit=Float64[], nsamples=20000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=Int[1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 21], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([50, 30, 500, 20], [50, 30, 500, 20]), maxtime=300.0, decayrate=1.0, totalsteps=20, verbose=false, hierarchical=(2, [8, 17], tuple()), method=(Tsit5(), true))
    rh = 50.0 .+ 10 * randn(ntrials)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials, hierarchical=(6, rh))
    tracesingle = Vector{Vector{Vector{Float64}}}(undef, length(R))
    df = Vector{DataFrame}(undef, length(R))
    k = 0
    for i in eachindex(R)
        tracesingle = [t[:, i] for t in trace]
        n = num_rates(transitions[i], R[i], S[i], insertstep[i]) + length(noisepriors[i])
        runit = rtarget[k+1:k+n]
        k += n
        df[i] = StochasticGene.make_traces_dataframe(tracesingle, interval, [runit; ones(length(runit)); repeat(runit, ntrials)], transitions[i], G[i], R[i], S[i], insertstep[i], prob_Gaussian, 4, "", true, true, tuple())
    end
    dfjoint = StochasticGene.make_traces_dataframe(trace, interval, [rtarget; ones(length(rtarget)); repeat(rtarget, ntrials)], transitions, G, R, S, insertstep, prob_Gaussian, 4, "", true, true, coupling)
    return df, dfjoint
end

function test_load(; traceinfo=(1.0, 1.0, -1, 1.0), datapath="data/inhibition/control/", label="", gene="MYC", datacond="gene", datatype=:trace, zscoretrace=true)
    data = load_data_trace(datapath, label, gene, datacond, traceinfo, datatype, 3, zscoretrace)
    return data
end


marg(p, dims) = dropdims(sum(p, dims=dims), dims=dims)


function extract_components(components)
    nR_values = []
    nG_values = []
    for c in components.modelcomponents
        push!(nR_values, c.nR)
        push!(nG_values, c.nG)
    end


    result = []
    for i in 1:length(nR_values)
        push!(result, nG_values[end-i+1])
        push!(result, nR_values[end-i+1])
    end

    return tuple(result...)
end

function test_a(Tc, GG)
    ac = StochasticGene.kolmogorov_forward(Tc', 1.0)
    ag = StochasticGene.kolmogorov_forward(GG', 1.0)
    return ac, ag
end


function test_distributions(r::Vector{Vector{type}}, coupling, transitions, G, R, S, insertstep, coupling_strength=[1.0, 1]) where {type<:Number}
    components = StochasticGene.make_components_TCoupled(coupling, transitions, G, R, S, insertstep)
    Tc = StochasticGene.make_mat_TC(components, r, coupling_strength)
    T, GG, Gt, Gs, IG, IR, IT = StochasticGene.make_matvec_C(components, r)
    ac = StochasticGene.kolmogorov_forward(Tc', 1.0)
    ag = StochasticGene.kolmogorov_forward(GG', 1.0)

    p0 = StochasticGene.normalized_nullspace(Tc)
    p01 = StochasticGene.normalized_nullspace(Tcr[1])
    p02 = StochasticGene.normalized_nullspace(Tcr[2])
    p0G = StochasticGene.normalized_nullspace(GG)

    p0R = reshape(p0, 2, 2, 2, 2)
    p01R = reshape(p01, 2, 2)
    p02R = reshape(p02, 2, 2, 2)
    pG = reshape(p0G, 2, 2)
    return p0R, p01R, p02R, pG
end

function recursive_sum(matrix, result, indices1, indices2, depth, max_depth)
    if depth == max_depth
        i, j, k, l = indices1
        i1, j1, k1, l1 = indices2
        idx1 = i + 2 * j + 4 * k + 8 * l + 1
        idx2 = i1 + 2 * j1 + 4 * k1 + 8 * l1 + 1
        result[i+2*k+1, i1+2*k1+1] += matrix[idx1, idx2]
    else
        for x in 0:1
            if depth < 4
                new_indices1 = copy(indices1)
                new_indices1[depth+1] = x
                recursive_sum(matrix, result, new_indices1, indices2, depth + 1, max_depth)
            else
                new_indices2 = copy(indices2)
                new_indices2[depth-4+1] = x
                recursive_sum(matrix, result, indices1, new_indices2, depth + 1, max_depth)
            end
        end
    end
end

function sum_overR2(matrix)
    result = zeros(4, 4)
    recursive_sum(matrix, result, [0, 0, 0, 0], [0, 0, 0, 0], 0, 8)
    return result / 4
end

function sum_overR(matrix)
    result = zeros(4, 4)

    for i in 0:1
        for j in 0:1
            for k in 0:1
                for l in 0:1
                    for i1 in 0:1
                        for j1 in 0:1
                            for k1 in 0:1
                                for l1 in 0:1
                                    result[i+2*k+1, i1+2*k1+1] += matrix[i+2j+4k+8l+1, i1+2j1+4k1+8l1+1]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    result / 4
end

function recursive_sum3(matrix, result, indices1, indices2, depth, max_depth, ranges)
    if depth == max_depth
        idx1 = sum(indices1[i] * prod(ranges[1:i-1]) for i in 1:length(indices1)) + 1
        idx2 = sum(indices2[i] * prod(ranges[1:i-1]) for i in 1:length(indices2)) + 1
        result[indices1[1]+2*indices1[3]+1, indices2[1]+2*indices2[3]+1] += matrix[idx1, idx2]
    else
        for x in 0:ranges[depth+1]-1
            if depth < length(indices1)
                indices1[depth+1] = x
                recursive_sum3(matrix, result, indices1, indices2, depth + 1, max_depth, ranges)
            else
                indices2[depth-length(indices1)+1] = x
                recursive_sum3(matrix, result, indices1, indices2, depth + 1, max_depth, ranges)
            end
        end
    end
end

function sum_overR3(matrix, ranges)
    half_length = div(length(ranges), 2)
    result_size = prod(ranges[1:half_length])
    result = zeros(eltype(matrix), result_size, result_size)
    indices1 = zeros(Int, half_length)
    indices2 = zeros(Int, half_length)
    recursive_sum3(matrix, result, indices1, indices2, 0, length(ranges), ranges)
    return result / result_size
end


function conditional_distributiona(joint_prob::Array, dist_index::Int)
    # Sum over the distribution index to get the marginal distribution
    marginal_prob = sum(joint_prob, dims=dist_index)

    # Ensure no zero probabilities to avoid division by zero
    marginal_prob[marginal_prob.==0] .= eps()

    # Expand the marginal to the shape of the joint distribution for broadcasting
    marginal_prob_expanded = repeat(marginal_prob, inner=(1, 1, 1, size(joint_prob, dist_index)))

    # Divide the joint distribution by the marginal to get the conditional distribution
    cond_prob = joint_prob ./ marginal_prob_expanded

    return cond_prob
end

function conditional_distribution(joint_prob::Array, dist_index::Int, cond_indices::Vector{Int})
    # Sum over the dimensions that are not part of the conditional variables
    dims_to_sum = setdiff(1:ndims(joint_prob), vcat(dist_index, cond_indices))
    marginal_prob = sum(joint_prob, dims=Tuple(dims_to_sum))  # Pass dims as a tuple

    # Ensure no zero probabilities to avoid division by zero
    marginal_prob[marginal_prob.==0] .= eps()

    # Expand the marginal to the shape of the joint distribution for broadcasting
    repeat_dims = ones(Int, ndims(joint_prob))
    for i in cond_indices
        repeat_dims[i] = size(joint_prob, i)
    end
    repeat_dims[dist_index] = size(joint_prob, dist_index)

    marginal_prob_expanded = repeat(marginal_prob, repeat_dims...)

    # Compute the conditional distribution by dividing the joint distribution by the expanded marginal
    cond_prob = joint_prob ./ marginal_prob_expanded

    return cond_prob
end



function split_matrix(mat::Matrix{Float64}, idx1::Vector{Int}, idx2::Vector{Int})
    return hcat([mat[idx1, j] for j in 1:size(mat, 2)]...,
        [mat[idx2, j] for j in 1:size(mat, 2)]...)
end

function split_matrix2(mat::Matrix{Float64})
    N, M = size(mat)
    half = div(N, 2)  # Divide N into two equal parts (assuming even N)

    return [mat[1:half, j] for j in 1:M],  # First row: first half
    [mat[half+1:N, j] for j in 1:M] # Second row: second half
end


function convert_to_m_dim(mat::Matrix{Float64})
    row1, row2 = split_matrix2(mat)
    return [[row1[j], row2[j]] for j in 1:size(mat, 2)]  # Each element is a vector of two vectors
end





"""
    compute_psis_loo(fits::Fit)

Compute PSIS-LOO (Pareto Smoothed Importance Sampling Leave-One-Out Cross-Validation).
Returns a tuple of (elpd_loo, elpd_loo_se) where:
- elpd_loo is the expected log pointwise predictive density for a new dataset
- elpd_loo_se is the standard error of elpd_loo

This implementation reuses the pointwise log-likelihoods stored in fits.lppd
that were computed during MCMC sampling for WAIC.
"""
function compute_psis_loo(fits::Fit)
    # Get pointwise log-likelihoods from fits.lppd
    n_samples = size(fits.param, 2)
    n_obs = length(fits.lppd)

    # Compute importance weights for each observation
    r_eff = compute_relative_eff(fits.param)  # Relative effective sample size

    # Compute LOO log-likelihoods using importance sampling
    loo_ll = zeros(n_obs)
    for i in 1:n_obs
        # Get pointwise log-likelihoods for this observation
        pointwise_ll = fits.lppd[i]

        # Compute importance weights for this observation
        weights = exp.(fits.param[i, :] .- maximum(fits.param[i, :]))
        weights .*= r_eff[i]
        weights ./= sum(weights)

        # Smooth weights using Pareto smoothing
        smoothed_weights = pareto_smooth_weights(weights)

        # Compute LOO estimate using smoothed weights
        loo_ll[i] = sum(smoothed_weights .* pointwise_ll)
    end

    # Compute elpd_loo and its standard error
    elpd_loo = sum(loo_ll)
    elpd_loo_se = sqrt(n_obs * var(loo_ll))

    return (elpd_loo, elpd_loo_se)
end

"""
    compute_relative_eff(params)

Compute relative effective sample size for importance sampling.
"""
function compute_relative_eff(params)
    n_samples = size(params, 2)
    ess = compute_ess(params)
    return ess ./ n_samples
end

"""
    pareto_smooth_weights(weights)

Apply Pareto smoothing to importance weights.
"""
function pareto_smooth_weights(weights)
    # Sort weights
    sorted_idx = sortperm(weights, rev=true)
    sorted_weights = weights[sorted_idx]

    # Find Pareto tail
    tail_idx = find_tail(sorted_weights)

    if tail_idx > 0
        # Fit Pareto distribution to tail
        k = fit_pareto(sorted_weights[tail_idx:end])

        # Replace tail with smoothed values
        smoothed = copy(weights)
        smoothed[sorted_idx[tail_idx:end]] = smooth_tail(sorted_weights[tail_idx:end], k)
        return smoothed
    end

    return weights
end

"""
    find_tail(weights)

Find the index where the Pareto tail begins.
"""
function find_tail(weights)
    n = length(weights)
    for i in 2:n
        if weights[i] / weights[1] < 0.1
            return i
        end
    end
    return 0
end

"""
    fit_pareto(weights)

Fit a Pareto distribution to the tail of weights.
Returns the shape parameter k.
"""
function fit_pareto(weights)
    n = length(weights)
    log_weights = log.(weights)
    k = 1 / (mean(log_weights) - log_weights[1])
    return k
end

"""
    smooth_tail(weights, k)

Smooth the tail of weights using the fitted Pareto distribution.
"""
function smooth_tail(weights, k)
    n = length(weights)
    smoothed = zeros(n)
    for i in 1:n
        smoothed[i] = weights[1] * (i / n)^(-1 / k)
    end
    return smoothed
end

"""
    read_measures_csv(filename::String)

Read the measures.csv file and return a DataFrame with the results.
"""
function read_measures_csv(filename::String)
    df = CSV.read(filename, DataFrame)
    # Sort by AIC to find best model
    sort!(df, :AIC)
    return df
end


function read_rates(rates_dir::String, model_name::String)
    if occursin("-h", model_name)
        rates_file = joinpath(rates_dir, "shared_trace-HBEC-$(model_name)_1.txt")
    else
        rates_file = joinpath(rates_dir, "rates_trace-HBEC-$(model_name)_1.txt")
    end
    readrates(rates_file)
end

function read_rates_params(rates_dir::String, model_name::String)
    rates_file = joinpath(rates_dir, "rates_trace-HBEC-$(model_name)_1.txt")
    param_file = joinpath(rates_dir, "param-stats_trace-HBEC-$(model_name)_1.txt")
    readrates(rates_file), readrates(param_file)
end

function decompose_nstate(model_name::String)
    parts = split(model_name, "_")
    G, R, S, insertstep = decompose_model(String(parts[end]))
    if G == 2
        transitions = ([1, 2], [2, 1])
    else
        transitions = ([1, 2], [2, 1], [2, 3], [3, 2])
    end
    return G, R, S, insertstep, transitions
end
"""
    evaluate_models_on_simulated_trace(measures_file::String, rates_dir::String)

Read measures.csv, find the best model by AIC, simulate a single long trace using only the shared parameters of the best model, then for all other models, read their rates, replace their noise/observation parameters with those from the simulation, and for hierarchical models use only the shared/hyper parameters for both rates and noise. Compute loglikelihood and AIC for all models on the simulated data using the existing loglikelihood functions. Returns a DataFrame with the new AIC and loglikelihood values for all models.
"""
function simulated_AIC(rates_dir::String; traceinterval=1.0, totaltime=1000.0, probfn=prob_Gaussian, n_noise=4, zeromedian=true, noisepriors=[0.0, 0.2, 1.0, 0.2], ntrials=40, nalleles=1)
    # using CSV, DataFrames

    # 1. Read measures file and sort by AIC
    measures_df = CSV.read(joinpath(rates_dir, "measures.csv"), DataFrame)
    sort!(measures_df, :AIC)

    # 2. Find the best model (lowest AIC)
    best_model = measures_df[1, :]
    best_model_name = String(best_model.Model)

    println(best_model_name)

    rates = read_rates(rates_dir, best_model_name)
    sim_noise = noisepriors
    if length(rates) >= n_noise
        rates[end-n_noise+1:end] .= sim_noise
    end

    G, R, S, insertstep, transitions = decompose_nstate(best_model_name)

    trace = simulate_trace_vector(rates, transitions, G, R, S, insertstep, traceinterval, totaltime, ntrials)

    nframes = round(Int, mean(size.(trace, 1)))  #mean number of frames of all traces

    data = TraceData{String,String,Tuple}("trace", "gene", traceinterval, (trace, 0.0, 0.0, nframes, 1.), Int[])

    results = DataFrame(Model=String[], AIC=Float64[], LogLikelihood=Float64[])
    for row in eachrow(measures_df)
        model_name = String(row.Model)

        rates = read_rates(rates_dir, model_name)

        # Replace noise parameters with those from the simulation
        if length(rates) >= n_noise
            rates[end-n_noise+1:end] .= sim_noise
        end

        G, R, S, insertstep, transitions = decompose_nstate(model_name)
        n = num_rates(transitions, R, S, insertstep)
        fittedparam = [collect(1:length(transitions)+2); collect(length(transitions)+R+1:n-1-max(0, S - 1))]

        model = load_model(data, rates, rates, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 0.1, Int[], 1., 0.01, prob_Gaussian, noisepriors, Tsit5(), tuple(), tuple(), nothing, zeromedian)

        param = transform_rates(rates[fittedparam], model)
        ll, _ = loglikelihood(param, data, model)
        n_params = length(fittedparam)
        aic = 2 * n_params - 2 * ll
        push!(results, (model_name, aic, ll))
    end
    sort!(results, :AIC), data
end

function dwelltime_AIC(rates_dir::String; bins=[collect(1:100), collect(0:100)], nhist=20, dttype=["ON", "OFF"], nalleles=1, onstates=[Int[], Int[]], total=10000000, tol=1e-6)
    # using CSV, DataFrames

    # 1. Read measures file and sort by AIC
    measures_df = CSV.read(joinpath(rates_dir, "measures.csv"), DataFrame)
    sort!(measures_df, :AIC)

    # 2. Find the best model (lowest AIC)
    best_model = measures_df[1, :]
    best_model_name = String(best_model.Model)

    println(best_model_name)

    G, R, S, insertstep, transitions = decompose_nstate(best_model_name)

    n = num_rates(transitions, R, S, insertstep)
    rates, _ = read_rates_params(rates_dir, best_model_name)

    dwelltimes = simulator(rates[1:n], transitions, G, R, S, insertstep, nalleles=nalleles, onstates=onstates[1], bins=bins[1], totalsteps=total, tol=tol)

    data = DwellTimeData("test", "test", bins, dwelltimes[2:end], dttype)

    results = DataFrame(Model=String[], AIC=Float64[], LogLikelihood=Float64[])
    for row in eachrow(measures_df)
        model_name = String(row.Model)

        rates, params = read_rates_params(rates_dir, model_name)

        n_params = length(params)
        println(n_params)

        G, R, S, insertstep, transitions = decompose_nstate(model_name)
        n = num_rates(transitions, R, S, insertstep)
        fittedparam = [collect(1:length(transitions)+2); collect(length(transitions)+R+1:n-1-max(0, S - 1))]

        model = load_model(data, rates, rates, fittedparam, tuple(), transitions, G, R, S, insertstep, "", nalleles, 0.1, onstates, 1., 0.1, prob_Gaussian, [], 1, tuple(), tuple(), nothing)

        param = transform_rates(rates[fittedparam], model)

        ll, _ = loglikelihood(param, data, model)

        aic = 2 * n_params - 2 * ll
        push!(results, (model_name, aic, ll))
    end
    sort!(results, :AIC), data
end

function dwelltime_matrices(rate_file::String, bins=[collect(1:100), collect(0:100)], onstates=[Int[], Int[]], dttype=["ON", "OFF"])
    model_name = remove_string(rate_file, "rates_trace-HBEC-", "_1.txt")
    rates = readrates(rate_file)
    data = DwellTimeData("test", "test", bins, ones(100), dttype)
    G, R, S, insertstep, transitions = decompose_nstate(model_name)
    println(G, R, S, insertstep)
    model = load_model(data, rates, rates, [], tuple(), transitions, G, R, S, insertstep, "", 1, 0.1, onstates, 1., 0.1, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    TI = make_mat(model.components.elementsTD[2], rates, model.components.nT)
    TA = make_mat(model.components.elementsTD[1], rates, model.components.nT)
    return TI, TA
end


# using StatsBase # For creating histograms later, if needed

"""
    process_reporter_times(df::DataFrame)

Analyzes a DataFrame to find columns named "Reporters..." and computes ON and OFF time durations.
`missing` values are ignored; the series is processed only up to the first `missing` value.
Typically, `missing` values are expected at the end of the series.

ON is defined as when the reporter value > 0. OFF is when reporter value == 0.
Edge effects are excluded:
- An OFF time is valid only if it's preceded by an ON state and followed by an ON state (ON -> OFF -> ON).
- An ON time is valid only if it's preceded by an OFF state and followed by an OFF state (OFF -> ON -> OFF).

# Arguments
- `df::DataFrame`: The input DataFrame. Reporter columns can contain numeric or missing values.

# Returns
- `Dict{String, Dict{String, Vector{Int}}}`: A dictionary where keys are reporter column names.
  Each value is another dictionary with keys "on_times" and "off_times",
  which are vectors of integers representing the durations of valid ON and OFF periods.
"""
function process_reporter_times(df::DataFrame)

    on_durations = Int[]
    off_durations = Int[]

    for rep_col_symb in reporter_col_symbols
        rep_col_name = string(rep_col_symb)
        series_with_potential_missing = df[!, rep_col_symb]

        # Determine the effective length of the series (up to the first 'missing')
        first_missing_index = findfirst(ismissing, series_with_potential_missing)

        local series::AbstractVector # Ensure series is defined in this scope
        if isnothing(first_missing_index)
            # No missing values, use the whole series
            series = series_with_potential_missing
        else
            # Truncate the series before the first missing value
            series = view(series_with_potential_missing, 1:(first_missing_index-1))
        end

        # 1. Determine states (True for ON, False for OFF) for the non-missing part
        # Now, 'val' is guaranteed not to be missing within the 'series'
        states = Bool[val > 0 for val in series]

        # 2. Identify blocks of consecutive states
        blocks = []
        if !isempty(states) # This check might be redundant given isempty(series) above, but safe
            current_block_start_idx = 1
            for i in 2:length(states)
                if states[i] != states[current_block_start_idx]
                    push!(blocks, (start_idx=current_block_start_idx, end_idx=i - 1, state=states[current_block_start_idx]))
                    current_block_start_idx = i
                end
            end
            push!(blocks, (start_idx=current_block_start_idx, end_idx=length(states), state=states[current_block_start_idx]))
        end

        num_blocks = length(blocks)

        # 3. Iterate through blocks to find valid ON/OFF durations
        if num_blocks >= 3
            for k in 2:(num_blocks-1)
                prev_block = blocks[k-1]
                current_block = blocks[k]
                next_block = blocks[k+1]

                duration = current_block.end_idx - current_block.start_idx + 1

                if current_block.state # Current block is ON
                    if !prev_block.state && !next_block.state # Preceded and followed by OFF
                        push!(on_durations, duration)
                    end
                else # Current block is OFF
                    if prev_block.state && next_block.state # Preceded and followed by ON
                        push!(off_durations, duration)
                    end
                end
            end
        end
        # results[rep_col_name] = Dict("on_times" => on_durations, "off_times" => off_durations)
    end
    return on_durations, off_durations, make_histogram(on_durations, normalize=true), make_histogram(off_durations, normalize=true)
end




# This is a conceptual example based on your hmm.jl and metropolis_hastings.jl
# You'll need to integrate these ideas into your actual codebase.

# Assumed necessary imports and struct definitions (HMMReporter, TComponents, etc.)
# using Distributions, LinearAlgebra, LogExpFunctions, CUDA, Distributed, StatsBase, LoopVectorization
# (Add other necessary structs and using statements from your project)

# --- Potentially Modified HMM functions (hmm.jl) ---

"""
set_b!(b_buffer::AbstractMatrix, trace_segment::AbstractVector, d::Vector{<:Distribution})

In-place version of set_b. Fills the pre-allocated b_buffer.
Assumes b_buffer has dimensions (N_states, length(trace_segment)).
"""
function set_b!(b_buffer::AbstractMatrix, trace_segment::AbstractVector, d::Vector{T_dist}) where {T_dist<:Distribution}
    N_states = length(d)
    T_len_segment = length(trace_segment)

    if size(b_buffer, 1) != N_states || size(b_buffer, 2) < T_len_segment
        # Or resize, or error more gracefully depending on strategy
        error("b_buffer dimensions are incompatible with N_states or trace_segment length.")
    end

    # Use a view if b_buffer is larger than needed for this specific trace_segment
    b_view = view(b_buffer, 1:N_states, 1:T_len_segment)

    for (t_idx, obs) in enumerate(trace_segment)
        for j_state in 1:N_states
            # Use @inbounds for slight performance gain if confident about bounds
            b_view[j_state, t_idx] = pdf(d[j_state], obs)
        end
    end
    # No return needed, b_buffer (via b_view) is modified
end

# Example for the multi-dimensional trace variant of set_b!
# function set_b!(b_buffer::AbstractMatrix, trace_segment_row::AbstractMatrix, d_vec_of_dist_vec::Vector{<:Vector{<:Distribution}})
#     N_states = length(d_vec_of_dist_vec[1]) # Assuming all inner vectors have same N_states
#     T_len_segment = size(trace_segment_row, 1) # Assuming trace_segment_row is T x num_features

#     if size(b_buffer, 1) != N_states || size(b_buffer, 2) < T_len_segment
#         error("b_buffer dimensions incompatible.")
#     end
#     b_view = view(b_buffer, 1:N_states, 1:T_len_segment)
#     fill!(b_view, 1.0) # Initialize with 1.0 for product

#     for t_idx in 1:T_len_segment
#         obs_features = view(trace_segment_row, t_idx, :) # Get the features for this time step
#         for j_state in 1:N_states
#             for i_feature in eachindex(d_vec_of_dist_vec) # Iterate over each feature's distributions
#                 b_view[j_state, t_idx] *= pdf(d_vec_of_dist_vec[i_feature][j_state], obs_features[i_feature])
#             end
#         end
#     end
# end


"""
forward!(alpha_buffer::AbstractMatrix, C_buffer::AbstractVector, a::AbstractMatrix, b_view::AbstractMatrix, p0::AbstractVector)

In-place version of the forward algorithm.
Fills pre-allocated alpha_buffer and C_buffer.
Assumes alpha_buffer has dimensions (N_states, T_len_segment)
Assumes C_buffer has length T_len_segment
Assumes b_view is already computed and has dimensions (N_states, T_len_segment)
"""
function forward!(alpha_buffer::AbstractMatrix, C_buffer::AbstractVector, a::AbstractMatrix, b_view::AbstractMatrix, p0::AbstractVector)
    N_states = size(a, 1)
    T_len_segment = size(b_view, 2)

    if size(alpha_buffer, 1) != N_states || size(alpha_buffer, 2) < T_len_segment || length(C_buffer) < T_len_segment
        error("Buffer dimensions are incompatible.")
    end

    # Use views for the buffers to only operate on the necessary part
    alpha_view = view(alpha_buffer, 1:N_states, 1:T_len_segment)
    C_view = view(C_buffer, 1:T_len_segment)

    # Initial step
    # Element-wise multiplication, then assign to the first column of alpha_view
    # Ensure p0 and b_view[:,1] are compatible for broadcasting if needed, or loop
    for j_state in 1:N_states
        alpha_view[j_state, 1] = p0[j_state] * b_view[j_state, 1]
    end

    current_sum = sum(view(alpha_view, :, 1))
    C_view[1] = 1.0 / max(current_sum, eps(Float64))
    # alpha_view[:, 1] .*= C_view[1] # In-place scaling
    for j_state in 1:N_states
        alpha_view[j_state, 1] *= C_view[1]
    end


    # Recursive step
    for t_idx in 2:T_len_segment
        fill!(view(alpha_view, :, t_idx), 0.0) # Zero out current time step in alpha_view
        for j_state in 1:N_states       # Current state
            sum_val = 0.0
            for i_state in 1:N_states   # Previous state
                sum_val += alpha_view[i_state, t_idx-1] * a[i_state, j_state]
            end
            alpha_view[j_state, t_idx] = sum_val * b_view[j_state, t_idx]
        end
        current_sum = sum(view(alpha_view, :, t_idx))
        C_view[t_idx] = 1.0 / max(current_sum, eps(Float64))
        # alpha_view[:, t_idx] .*= C_view[t_idx] # In-place scaling
        for j_state in 1:N_states
            alpha_view[j_state, t_idx] *= C_view[t_idx]
        end
    end
    # No explicit return needed, buffers are modified.
    # However, returning the sum of log(C) is often useful for log-likelihood.
    # This function would be part of a larger log-likelihood calculation.
end


# --- Modified _ll_hmm (conceptual) ---
# This is for the hierarchical case where `a` and `p0` are shared,
# but `noiseindividual` means `d` (and thus `b`) is trace-specific.

# Helper function to get max trace length if not already available
function get_max_trace_length(traces::Vector{<:AbstractArray})
    isempty(traces) && return 0
    # Assuming traces[i] is T_frames x N_features or just T_frames
    return maximum(size(t, 1) for t in traces)
end

"""
_ll_hmm_preallocated(
    a::AbstractMatrix, p0::AbstractVector,
    noiseindividual::Vector{<:Vector{Float64}}, # Vector of noise params for each trace
    reporter::HMMReporter, # Or Vector{HMMReporter} if it varies per trace (unlikely)
    traces::Vector{<:AbstractArray}, # Vector of trace matrices/vectors
    # Pre-allocated buffers:
    b_buffer::AbstractMatrix,
    alpha_buffer::AbstractMatrix,
    C_buffer::AbstractVector,
    d_buffer::Vector{<:Distribution} # Buffer for distributions if `set_d` allocates
)

Calculates log-likelihood using pre-allocated buffers.
"""
function _ll_hmm_preallocated(
    a::AbstractMatrix, p0::AbstractVector,
    noiseindividual::Vector{<:Vector{Float64}},
    reporter::HMMReporter, # Assuming a single reporter structure for all traces for simplicity
    traces::Vector{<:AbstractArray},
    b_buffer::AbstractMatrix,
    alpha_buffer::AbstractMatrix,
    C_buffer::AbstractVector,
    d_buffer::Vector{T_dist} # Buffer for the distributions
) where {T_dist<:Distribution}

    total_log_likelihood = 0.0
    # logpredictions = Vector{Float64}(undef, length(traces)) # If you need individual logpreds

    N_states = size(a, 1)

    for i in eachindex(traces)
        current_trace = traces[i]
        current_trace_len = size(current_trace, 1) # Assuming trace is T x features or T-vector

        # 1. Set distributions `d` for the current trace's noise parameters
        # If set_d itself allocates, you might need a set_d!
        # For now, assume set_d is relatively cheap or d_buffer is used by it.
        # This example assumes set_d! populates d_buffer
        # set_d!(d_buffer, noiseindividual[i], reporter)
        # If set_d returns a new vector of distributions:
        current_d = set_d(noiseindividual[i], reporter) # Original call from your code

        # 2. Compute emission probabilities `b` into b_buffer
        # Pass only the relevant part of the trace if it's a matrix (e.g. trace[i] for vector of traces)
        set_b!(b_buffer, current_trace, current_d)
        b_view_for_current_trace = view(b_buffer, 1:N_states, 1:current_trace_len)

        # 3. Run forward algorithm using alpha_buffer and C_buffer
        forward!(alpha_buffer, C_buffer, a, b_view_for_current_trace, p0)

        # Calculate log-likelihood for this trace
        # C_view_for_current_trace = view(C_buffer, 1:current_trace_len) # Not strictly needed if sum is done carefully
        current_trace_log_likelihood = -sum(log.(view(C_buffer, 1:current_trace_len)))
        total_log_likelihood += current_trace_log_likelihood
        # logpredictions[i] = current_trace_log_likelihood
    end

    return total_log_likelihood #, logpredictions (if needed)
end


# --- How to integrate into your metropolis_hastings or main MCMC loop ---

# Inside the function that calls _ll_hmm (e.g., within metropolis_hastings or a worker function):
# This would be done ONCE per worker process, before the MCMC loop starts.

# Example:
# function run_single_mcmc_chain(data, model, options, trace_data_for_this_chain)
# ... (setup a, p0, initial params etc.) ...

# Determine max_trace_len and N_states from data/model
# max_T = get_max_trace_length(trace_data_for_this_chain) # Or data.trace[1] if all traces are there
# N_states = model.components.nT # Or derive appropriately

# Pre-allocate buffers based on max possible dimensions
# b_buffer = Matrix{Float64}(undef, N_states, max_T)
# alpha_buffer = Matrix{Float64}(undef, N_states, max_T)
# C_buffer = Vector{Float64}(undef, max_T)
# d_buffer = Vector{Distribution{Univariate, Continuous}}(undef, N_states) # Example

# Inside your MCMC loop (e.g., in mhstep, which calls loglikelihood, which calls _ll_hmm):
# ...
# current_loglik, _ = _ll_hmm_preallocated(
#                             a, p0,
#                             current_noise_individual_params, # Extracted for current MCMC sample
#                             model.reporter,
#                             trace_data_for_this_chain, # Or data.trace[1]
#                             b_buffer, alpha_buffer, C_buffer, d_buffer
#                         )
# ...
# end


#= Feature/regression helper moved to `test_features.jl`.

"""
    test_run_spec_roundtrip(; tmpdir=mktempdir())

Round-trip test for the JLD2 run_spec serialization.  Builds a run_spec dict that mirrors
the exact types produced by fit(...) for a coupled tracejoint model, writes it to a JLD2
companion file via write_run_spec_jld2, reads it back via read_run_spec, then asserts every
field is identical (===, ==, or ≈ as appropriate).

Usage:
    StochasticGene.test_run_spec_roundtrip()
"""
function test_run_spec_roundtrip(; tmpdir=mktempdir())
    run_spec = Dict{Symbol, Any}(
        :key            => "11",
        :datapath       => "data/3Prime_gene_enhancer/including_background/short",
        :TransitionType => "nstate",
        :maxtime        => 30.0,
        :inlabel        => "tracejoint-HBEC-nstate_enhancer-gene11",
        :label          => "tracejoint-HBEC-nstate_enhancer-gene11",
        :nchains        => 1,
        :samplesteps    => 100000,
        :datatype       => "tracejoint",
        :warmupsteps    => 0,
        :cell           => "HBEC",
        :resultfolder   => "3Prime-coupled-test",
        :propcv         => 0.05,
        :ratetype       => "median",
        :infolder       => "3Prime-coupled-2026-03-11",
        :datacol        => 3,
        :gene           => "MYC",
        :root           => ".",
        :decayrate      => 1.0,
        :transitions    => (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])),
        :G              => (3, 3),
        :R              => (3, 3),
        :S              => (0, 0),
        :insertstep     => (1, 1),
        :coupling       => ((1, 2), [(1, 1, 2, 1)], [:inhibit]),
        :priormean      => [0.001, 0.001, 0.01, 0.01, 0.1, 0.15, 0.15, 0.15, 0.03165055618995184, 0.0, 0.1, 0.5, 0.15, 0.01, 0.01, 0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 0.0, 0.1, 0.9, 0.2, -0.1],
        :priorcv        => [1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.01, 0.1, 10.0],
        :fittedparam    => [14, 27],
        :fixedeffects   => (),
        :noisepriors    => ([0.0, 0.1, 0.5, 0.15], [0.0, 0.1, 0.9, 0.2]),
        :hierarchical   => (),
        :datacond       => ["enhancer", "gene"],
        :elongationtime => (20.0, 5.0),
        :traceinfo      => (1.6666666666666667, 1.0, -1),
        :zeromedian     => Bool[1, 1],
        :probfn         => (prob_Gaussian, prob_Gaussian),
        :method         => Tsit5(),
        :prerun         => 0.0,
    )

    # Write and read back via JLD2
    fake_toml = joinpath(tmpdir, "info_roundtrip.toml")
    write_run_spec_jld2(fake_toml, run_spec)
    loaded = read_run_spec(fake_toml)

    failures = String[]

    for k in keys(run_spec)
        orig  = run_spec[k]
        back  = get(loaded, k, missing)
        if back === missing
            push!(failures, "  MISSING key :$k")
            continue
        end
        ok = try
            if orig isa AbstractFloat || (orig isa AbstractArray && eltype(orig) <: AbstractFloat)
                orig ≈ back
            else
                orig == back
            end
        catch
            false
        end
        ok || push!(failures, "  MISMATCH :$k\n    orig : $(repr(orig))\n    back : $(repr(back))")
    end

    if isempty(failures)
        println("run_spec round-trip: ALL $(length(run_spec)) fields OK  (tmpdir=$tmpdir)")
    else
        println("run_spec round-trip: $(length(failures)) FAILURE(S):")
        foreach(println, failures)
    end
    return isempty(failures)
end
=#

"""
    profile_trace_prediction(info_jld2::String, rates_file::String; ratetype="median")

Profile a single `make_traces_dataframe_key` call to find allocation hotspots.
Runs the call once to warm up JIT, then profiles allocations on the second call.

Reports:
- Wall time and allocations for the full call (@time)
- Top allocation sites via Profile.Allocs (requires Julia 1.8+)

Usage:
    StochasticGene.profile_trace_prediction(
        "results/5Prime-coupled-2026-03-11/info_2435.jld2",
        "results/5Prime-coupled-2026-03-11/rates_2435.txt")
"""
function profile_trace_prediction(info_jld2::String, rates_file::String; ratetype::String="median")
    info = read_run_spec(info_jld2)
    r    = readrates(rates_file, get_row(ratetype))
    traceinfo = info[:traceinfo]
    data = load_data_trace(info[:datapath], "", "", info[:datacond],
        (Float64(traceinfo[1]), traceinfo[2], traceinfo[3], 1.0, 0.0),
        :trace, info[:datacol], info[:zeromedian])

    println("=== Warmup (JIT) ===")
    make_traces_dataframe_key(data, r, info)

    println("\n=== Timed run ===")
    @time make_traces_dataframe_key(data, r, info)

    bytes = @allocated make_traces_dataframe_key(data, r, info)
    println("\n=== Allocations (timed run, post-JIT) ===")
    println("  Total: ", round(bytes / 1024^3, digits=2), " GiB")
    println("\nTo get per-site breakdown, run at the REPL:")
    println("  using Profile")
    println("  Profile.Allocs.@profile sample_rate=1 StochasticGene.make_traces_dataframe_key(data, r, info)")
    println("  PProf.Allocs.pprof()  # requires PProf.jl")
    return nothing
end

#= Legacy benchmark implementations moved to `benchmarks.jl`.

# ════════════════════════════════════════════════════════════════════════════════════
# SECTION 5: DEVELOPER BENCHMARKS & PERFORMANCE ANALYSIS
# ════════════════════════════════════════════════════════════════════════════════════
# Performance testing infrastructure for inference, gradients, and algorithm comparison.
# NOT called by Pkg.test(); designed for development, profiling, and performance validation.
# Useful for benchmarking AD methods (Zygote, ForwardDiff, finite-diff) and inference strategies.

# --- Inference benchmark (MH / NUTS / ADVI wall-clock comparison) ---
# Interactive: `using StochasticGene; benchmark_inference_setup_parallel_workers(8)` then
# `scen = benchmark_inference_simrna_small()` and `benchmark_inference_run_mh(scen; maxtime=600.0)`.
# Or start Julia with extra workers: `julia --project=. -p 8`.

"""
    benchmark_inference_simrna_small(; kwargs...) -> NamedTuple

Build a small **GM + RNA histogram** [`RNAData`](@ref) and [`GMmodel`](@ref) for comparing
[`run_mh`](@ref), [`run_nuts_fit`](@ref), and [`run_advi`](@ref) on the same synthetic problem.

# Returns
`(; data, model, meta)` — pass `scenario` to [`benchmark_inference_run_mh`](@ref), etc.

# Keyword arguments
- `rtarget`, `transitions`, `G`, `nRNA`, `nalleles`, `fittedparam`, `fixedeffects`, `rinit`, `totalsteps`: simulation / model layout (defaults match [`test_fit_simrna`](@ref)-style smoke tests).
- `seed`: passed to `Random.seed!` when not `nothing` (default `42`).

# Example
```julia
using StochasticGene
scen = benchmark_inference_simrna_small(seed=1)
benchmark_inference_run_mh(scen; nchains=1, maxtime=30.0)  # laptop smoke
```
"""
function benchmark_inference_simrna_small(;
    rtarget=[0.33, 0.19, 2.5, 1.0],
    transitions=([1, 2], [2, 1]),
    G::Int=2,
    nRNA::Int=100,
    nalleles::Int=2,
    fittedparam=[1, 2, 3],
    fixedeffects=tuple(),
    rinit=[0.1, 0.1, 0.1, 1.0],
    totalsteps::Int=100_000,
    seed::Union{Nothing,Int}=42,
)
    seed !== nothing && Random.seed!(seed)
    h = simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=totalsteps, nalleles=nalleles)[1]
    data = RNAData{typeof(nRNA),typeof(h)}("", "", nRNA, h, 1.0, [1])
    model = load_model(
        data, rinit,
        prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0),
        fittedparam, fixedeffects, transitions, G, 0, 0, 0, "", nalleles, 10.0, Int[], rtarget[end], 0.02,
        prob_Gaussian, [], 1, tuple(), tuple(), nothing,
    )
    meta = (; rtarget, transitions, G, nRNA, nalleles, fittedparam, fixedeffects, rinit, totalsteps, seed)
    return (; data, model, meta)
end

"""
    benchmark_inference_trace_gr2r2(; kwargs...) -> NamedTuple

Simulated **trace** data and a **GRSM** model with **G = 2**, **R = 2**, **S = 0**, **insertstep = 1**
(same layout as [`test_fit_trace`](@ref)). Use with [`benchmark_inference_run_mh`](@ref),
[`benchmark_inference_run_nuts_parallel`](@ref), [`benchmark_inference_run_advi`](@ref).

For NUTS on traces, pass **`steady_state_solver=:augmented`** to [`benchmark_inference_run_nuts_parallel`](@ref)
(the RNA-histogram default `:default` is for steady-state RNA likelihoods).

# Returns
`(; data, model, meta)` with `TraceData` and a loaded `GRSMmodel`.

# Keyword arguments
- `traceinfo`, `transitions`, `rtarget`, `totaltime`, `ntrials`, `fittedparam`, `propcv`, `noisepriors`, `zeromedian`, `initprior`: same roles as [`test_fit_trace`](@ref). Default **`propcv=0.01`** matches [`fit`](@ref); increase only if MH acceptance is too high.
- `seed`: passed to `Random.seed!` when not `nothing` (default `42`).

# Example
```julia
using StochasticGene
scen = benchmark_inference_trace_gr2r2(seed=1, totaltime=500.0, ntrials=4)
benchmark_inference_run_mh(scen; nchains=1, maxtime=60.0)
benchmark_inference_run_nuts_parallel(scen; nchains=2, n_samples=200, steady_state_solver=:augmented)
```
"""
function benchmark_inference_trace_gr2r2(;
    traceinfo=(1.0, 1.0, -1, 1.0, 0.5),
    G::Int=2,
    R::Int=2,
    S::Int=0,
    insertstep::Int=1,
    transitions=([1, 2], [2, 1]),
    rtarget=[0.01, 0.01, 0.1, 0.15, 0.1, 1.0, 0, 0.1, 1., 0.1],
    totaltime::Float64=4000.0,
    ntrials::Int=10,
    fittedparam=[1,2,3,4,5,7],
    propcv::Float64=0.01,
    noisepriors=[0.0, 0.1, 1.0, 0.1],
    zeromedian::Bool=true,
    initprior::Float64=0.1,
    seed::Union{Nothing,Int}=42,
)
    seed !== nothing && Random.seed!(seed)
    tracer = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, traceinfo[1], totaltime, ntrials)
    trace, tracescale = zero_median(tracer, zeromedian)
    nframes = round(Int, mean(size.(trace, 1)))
    if length(traceinfo) > 3 && traceinfo[4] != 1.0
        weight = set_trace_weight(traceinfo)
        background = set_trace_background(traceinfo)
    else
        weight = 0.0
        background = 0.0
    end
    data = TraceData{String,String,Tuple}("trace", "gene", traceinfo[1], (trace, background, weight, nframes, tracescale), Int[])
    elongationtime = mean_elongationtime(rtarget, transitions, R)
    priormean = [fill(0.01, length(transitions)); initprior; fill(R / elongationtime, R); fill(0.05, max(0, S - insertstep + 1)); 1.0; noisepriors]
    priorcv = [fill(1.0, length(transitions)); 0.1; fill(0.1, R); fill(0.1, max(0, S - insertstep + 1)); 1.0; [0.5, 0.5, 0.1, 0.1]]
    rinit = !isempty(rtarget) ? set_rinit(rtarget, priormean) : set_rinit(rtarget, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), tuple(), nothing)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, priorcv, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, Tsit5(), tuple(), tuple(), nothing, zeromedian)
    meta = (; traceinfo, G, R, S, insertstep, transitions, rtarget, totaltime, ntrials, fittedparam, propcv, zeromedian, seed)
    return (; data, model, meta)
end

"""
    benchmark_inference_trace_coupled_3x3(; kwargs...) -> NamedTuple

Simulated **coupled joint-trace** data and a **GRSM** model with **two units**, each
**`(G,R,S,insertstep) = (3,3,0,1)`**, 3-state promoter transitions per unit, and reciprocal coupling
`(1,2)` with one hidden connection — the same geometry as [`test_benchmark_trace_joint_fit_stacks`](@ref)
defaults. Use with [`benchmark_inference_run_mh`](@ref), [`benchmark_inference_run_nuts_parallel`](@ref),
[`benchmark_trace_finitediff_gradient`](@ref), [`benchmark_trace_zygote_subset_gradient`](@ref), etc. Do **not** use
[`benchmark_trace_zygote_gradient`](@ref) (full-parameter Zygote) on this scenario unless you pass **`allow_coupled_zygote=true`**
— it can OOM. Use [`hmm_checkpoint_steps`](@ref) and subset-gradient benchmarks for coupled traces.

# Returns
`(; data, model, meta)` with `TraceData` (`"tracejoint"`) and [`load_model`](@ref) `GRSMmodel`.
Default **`coupled_stack=:full`** (`TCoupledFullComponents`); pass **`coupled_stack=:legacy`** to match the
RG / Kronecker stack used in [`test_benchmark_trace_joint_fit_stacks`](@ref) timing comparisons.

# Keyword arguments
- `seed`: `Random.seed!` when not `nothing` (default `42`).
- `coupling`, `G`, `R`, `S`, `insertstep`, `transitions`, `rtarget`: coupled layout (defaults match
  [`test_benchmark_trace_joint_fit_stacks`](@ref); `rtarget` length is `num_rates(...) + ncoupling`).
- `totaltime`, `ntrials`, `interval`, `trace_specs`, `noisepriors`: simulation and [`TraceData`](@ref) construction
  (default `trace_specs` is [`default_trace_specs_for_coupled`](@ref)).
- `fittedparam`, `propcv`, `method`, `zeromedian`, `coupled_stack`: passed to [`load_model`](@ref).

# Example
```julia
using StochasticGene
scen = benchmark_inference_trace_coupled_3x3(seed=1, totaltime=500.0, ntrials=2)
benchmark_inference_run_mh(scen; nchains=1, maxtime=60.0)
```
"""
function benchmark_inference_trace_coupled_3x3(;
    seed::Union{Nothing,Int}=42,
    coupling=((1, 2), [(1, 2, 2, 1)], [:free]),
    G=(3, 3), R=(3, 3), S=(0, 0), insertstep=(1, 1),
    transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])),
    rtarget::Vector{Float64}=Float64[
        0.03, 0.1, 0.5, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2,
        0.0, 0.1, 0.5, 0.15,
        0.03, 0.1, 0.5, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2,
        0.0, 0.1, 0.9, 0.2,
        -0.4,
    ],
    totaltime::Float64=800.0,
    ntrials::Int=4,
    interval::Float64=1.0,
    noisepriors=([0., .1, 1., .1], [0., .1, 1., .1]),
    trace_specs=nothing,
    fittedparam=Int[19],
    propcv::Float64=0.2,
    method=Tsit5(),
    coupled_stack::Symbol=:full,
    zeromedian::Bool=true,
)
    seed !== nothing && Random.seed!(seed)
    trace_specs_eff = trace_specs === nothing ? StochasticGene.default_trace_specs_for_coupled((interval, 1.0, -1.0), [false, false], 2) : trace_specs
    units = [spec.unit for spec in trace_specs_eff]
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling,
                                  interval, totaltime, ntrials; noiseparams=[4, 4])
    data = TraceData("tracejoint", "test", interval, (trace, [], fill(0.0, length(units)), 1), units)
    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors,
                              mean_elongationtime(rtarget, transitions, R), tuple(), coupling, nothing)
    rinit = rtarget
    nr = num_rates(transitions, R, S, insertstep)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep,
                       "", 1, 10.0, Int[], rtarget[nr], propcv, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing,
                       zeromedian; coupled_stack=coupled_stack)
    meta = (; coupling, G, R, S, insertstep, transitions, rtarget, totaltime, ntrials, interval,
            fittedparam, propcv, coupled_stack, zeromedian, seed)
    return (; data, model, meta)
end

"""
    benchmark_scenario_coupled_3x3_10traces_220frames(; kwargs...) -> NamedTuple

Convenience wrapper around [`benchmark_inference_trace_coupled_3x3`](@ref) with **`ntrials=10`**,
**`totaltime=220`**, **`interval=1.0`**, targeting on the order of **220 frames per trace** (exact length follows the simulator).

Additional `kwargs` are forwarded (e.g. `seed`, `coupled_stack`, `fittedparam`, `rtarget`).
"""
function benchmark_scenario_coupled_3x3_10traces_220frames(; kwargs...)
    benchmark_inference_trace_coupled_3x3(; ntrials=10, totaltime=220.0, interval=1.0, kwargs...)
end

"""
    benchmark_inference_trace_coupled_3x3_g3r0(; kwargs...) -> NamedTuple

Simulated **three-unit coupled joint-trace** scenario: **two** observed **(G,R,S,insertstep) = (3,3,0,1)**
reporter units plus a **hidden** **(G,R,S) = (3,0,0)** telegraph (**R = 0** on that unit). Same layout as
[`test_fit_tracejoint_3unit`](@ref) but with **G = (3,3,3)** and **R = (3,3,0)** instead of the smaller 2+1 reporter geometry.

The third unit uses **`insertstep = 1`** (matching [`test_fit_tracejoint_3unit`](@ref)); **`insertstep = 0`** with **`R = 0`**
mis-sizes the auto-generated prior vector relative to [`num_rates`](@ref).

- **Coupling** (default): `unit_model = (1,2,3)`, connections `(3,1,1,2)` and `(3,3,2,2)`, **`[:inhibit,:inhibit]`**
  (hidden unit 3 affects units 1 and 2).
- **Traces**: only units **1** and **2** (`data.units`; unit 3 is latent). Uses [`default_trace_specs_for_coupled`](@ref)
  with [`n_observed_trace_units`](@ref)`(coupling) == 2` when `trace_specs === nothing`.
- **Parameters**: default **`fixedeffects`** ties all unit-3 rates to one master and both coupling strengths to one
  (same pattern as [`test_fit_tracejoint_3unit`](@ref)); override with **`fixedeffects`** / **`fittedparam`** /
  **`rtarget`** for other experiments.

As with [`benchmark_inference_trace_coupled_3x3`](@ref), [`benchmark_trace_zygote_gradient`](@ref) is **opt-in only**
(`allow_coupled_zygote=true`); use **finite differences** for AD benchmarks on coupled traces ([`benchmark_trace_finitediff_gradient`](@ref)).

# Returns
`(; data, model, meta)` with `TraceData` (`"tracejoint"`) and [`load_model`](@ref) `GRSMmodel`.

# Example
```julia
using StochasticGene
scen = benchmark_inference_trace_coupled_3x3_g3r0(seed=1, totaltime=500.0, ntrials=2)
benchmark_inference_run_mh(scen; nchains=1, maxtime=120.0)
```
"""
function benchmark_inference_trace_coupled_3x3_g3r0(;
    seed::Union{Nothing,Int}=42,
    coupling=((1, 2, 3), [(3, 1, 1, 2), (3, 3, 2, 2)], [:inhibit, :inhibit]),
    G=(3, 3, 3),
    R=(3, 3, 0),
    S=(0, 0, 0),
    insertstep=(1, 1, 1),
    transitions=(
        ([1, 2], [2, 1], [2, 3], [3, 2]),
        ([1, 2], [2, 1], [2, 3], [3, 2]),
        ([1, 2], [2, 1], [2, 3], [3, 2], [1, 3], [3, 1]),
    ),
    rtarget::Union{Nothing,Vector{Float64}}=nothing,
    units=[1, 2],
    totaltime::Float64=800.0,
    ntrials::Int=4,
    interval::Float64=1.0,
    noisepriors=([0., .1, 1., .1], [0., .1, 1., .1], Float64[]),
    trace_specs=nothing,
    fixedeffects=nothing,
    fittedparam=nothing,
    propcv::Float64=0.2,
    method=Tsit5(),
    coupled_stack::Symbol=:full,
    zeromedian::Bool=true,
)
    seed !== nothing && Random.seed!(seed)
    nrates_per_unit = num_rates(transitions, R, S, insertstep)
    ncpl = ncoupling(coupling)
    n_noise = sum(length.(noisepriors))
    n_total = sum(nrates_per_unit) + n_noise + ncpl
    unit3_rate_start = sum(nrates_per_unit[1:2]) + sum(length.(noisepriors[1:2])) + 1
    unit3_rate_end = unit3_rate_start + nrates_per_unit[3] - 1
    coupling_start = sum(nrates_per_unit) + n_noise + 1
    if rtarget === nothing
        u1r = Float64[0.03, 0.1, 0.5, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2]
        u1n = Float64[0.0, 0.1, 0.5, 0.15]
        u2r = Float64[0.03, 0.1, 0.5, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2]
        u2n = Float64[0.0, 0.1, 0.9, 0.2]
        u3r = Float64[0.12, 0.08, 0.15, 0.1, 0.1, 0.1, 0.1, 0.2]
        rtarget = vcat(u1r, u1n, u2r, u2n, u3r, fill(-0.2, ncpl))
    else
        length(rtarget) == n_total || throw(ArgumentError("rtarget must have length $n_total (sum(num_rates)+noise+coupling)"))
    end
    if fixedeffects === nothing
        fixedeffects = (collect(unit3_rate_start:unit3_rate_end), collect(coupling_start:coupling_start + ncpl - 1))
    end
    if fittedparam === nothing
        fittedparam = Int[unit3_rate_start, coupling_start]
    end
    trace_specs_eff = if trace_specs === nothing
        zm = zeromedian isa Bool ? fill(zeromedian, length(units)) : zeromedian
        StochasticGene.default_trace_specs_for_coupled((interval, 1.0, -1.0), zm, units)
    else
        trace_specs
    end
    units_eff = [spec.unit for spec in trace_specs_eff]
    length(units_eff) == length(units) || throw(ArgumentError("trace_specs units must match observed unit count"))
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling,
                                  interval, totaltime, ntrials; observed_units=units_eff, noiseparams=[4, 4, 0])
    data = TraceData("tracejoint", "test", interval, (trace, [], fill(0.0, length(units_eff)), 1), units_eff)
    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors,
                              mean_elongationtime(rtarget, transitions, R), tuple(), coupling, nothing)
    rinit = rtarget
    block_ends = cumsum([nrates_per_unit[i] + length(noisepriors[i]) for i in eachindex(R)])
    block_starts = vcat(1, block_ends[1:end-1] .+ 1)
    decayrate = tuple((rtarget[block_starts[i] + nrates_per_unit[i] - 1] for i in eachindex(R))...)
    model = load_model(data, rinit, priormean, fittedparam, fixedeffects, transitions, G, R, S, insertstep,
                       "", 1, 10.0, Int[], decayrate, propcv, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing,
                       zeromedian; coupled_stack=coupled_stack)
    meta = (; coupling, G, R, S, insertstep, transitions, rtarget, totaltime, ntrials, interval, units=units_eff,
            fittedparam, fixedeffects, propcv, coupled_stack, zeromedian, seed)
    return (; data, model, meta)
end

"""
    benchmark_inference_ensure_workers(n::Int)

Require at least `n` Julia **worker** processes for parallel chains (`run_mh(..., nchains)` and
[`benchmark_inference_run_nuts_parallel`](@ref)). Does not add workers.

Throws `ArgumentError` if `nworkers() < n`. Start Julia with e.g. `julia -p 8`, or run
[`benchmark_inference_setup_parallel_workers`](@ref) once per session.
"""
function benchmark_inference_ensure_workers(n::Int)
    n < 1 && throw(ArgumentError("n must be >= 1"))
    nworkers() >= n && return nothing
    throw(ArgumentError(
        "Need at least $n worker processes (nworkers()=$(nworkers())). " *
        "Start with e.g. `julia -p $n` or run `using Distributed; addprocs($n)` then " *
        "`@everywhere using StochasticGene`, or call `benchmark_inference_setup_parallel_workers($n)`.",
    ))
end

"""
    benchmark_inference_setup_parallel_workers(n::Int) -> Int

Add workers until `nworkers() >= n` and load `StochasticGene` on all workers (`@everywhere using`).
Returns `nworkers()` after setup. Safe for interactive benchmarking on a laptop; optional for CI.
"""
function benchmark_inference_setup_parallel_workers(n::Int)
    n < 1 && throw(ArgumentError("n must be >= 1"))
    if nworkers() < n
        addprocs(n - nworkers())
    end
    # `@everywhere` cannot run inside `include`/`test.jl`; load on each worker explicitly.
    for w in workers()
        remotecall_eval(Main, w, :(using StochasticGene))
    end
    return nworkers()
end

function _benchmark_ess_summary(measures::Measures, elapsed_sec::Float64)
    ess = measures.ess
    min_ess = minimum(ess)
    mean_ess = mean(ess)
    max_rhat = maximum(measures.rhat)
    return (;
        min_ess=min_ess,
        mean_ess=mean_ess,
        max_rhat=max_rhat,
        ess_per_sec_min=min_ess / elapsed_sec,
        ess_per_sec_mean=mean_ess / elapsed_sec,
    )
end

"""
    benchmark_inference_run_mh(data, model; kwargs...)
    benchmark_inference_run_mh(scenario; kwargs...)

Parallel Metropolis–Hastings with [`run_mh`](@ref). Requires [`benchmark_inference_ensure_workers`](@ref)(`nchains`) unless `nchains==1`.

# Keywords
- `nchains`, `samplesteps`, `warmupsteps`, `maxtime`: passed to [`MHOptions`](@ref) (`maxtime` in seconds, total wall budget per chain for warmup+sampling).
- `mh_seed`: optional `Random.seed!` before running.

Proposal scaling is **not** set here — it lives on `model` as **`propcv`** (see [`load_model`](@ref)). If you see “Low acceptance” / huge `max|Δparam|` / `prior=-Inf` from [`metropolis_hastings`](@ref), **rebuild** the scenario with a smaller `propcv` (e.g. [`benchmark_inference_trace_gr2r2`](@ref)(`propcv=0.005`) for trace GRSM). [`fit`](@ref) defaults to `propcv=0.01`.

# Returns
`NamedTuple` with `elapsed_sec`, `fits`, `stats`, `measures`, and ESS summary fields (`min_ess`, `ess_per_sec_min`, …).
"""
function benchmark_inference_run_mh(
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel;
    nchains::Int=1,
    samplesteps::Int=10_000,
    warmupsteps::Int=5_000,
    maxtime::Float64=600.0,
    mh_seed::Union{Nothing,Int}=nothing,
)
    nchains > 1 && benchmark_inference_ensure_workers(nchains)
    mh_seed !== nothing && Random.seed!(mh_seed)
    opts = MHOptions(samplesteps, warmupsteps, 0, maxtime_seconds(maxtime), 1.0, 2.0)
    elapsed_sec = @elapsed begin
        fits, stats, measures = run_mh(data, model, opts, nchains)
    end
    ess = _benchmark_ess_summary(measures, elapsed_sec)
    return (; elapsed_sec, fits, stats, measures, ess...)
end

function benchmark_inference_run_mh(scenario::NamedTuple; kwargs...)
    hasproperty(scenario, :data) && hasproperty(scenario, :model) ||
        throw(ArgumentError("scenario must be a NamedTuple with :data and :model (e.g. from benchmark_inference_simrna_small)"))
    benchmark_inference_run_mh(scenario.data, scenario.model; kwargs...)
end

"""
    benchmark_inference_run_nuts_parallel(data, model; kwargs...)
    benchmark_inference_run_nuts_parallel(scenario; kwargs...)

For **`nchains > 1`**, runs independent [`run_nuts_fit`](@ref) chains in parallel (`Distributed.@spawn`), then
[`merge_fit`](@ref) and pooled diagnostics. For **`nchains == 1`**, runs **in the main process** (no `@spawn`) so **ForwardDiff**, progress bars, and verbose output behave reliably (workers can break or wrap errors).

Use **`steady_state_solver=:default`** to align log-likelihood with MH for RNA histogram benchmarks (see [`test_compare_mh_nuts_posterior`](@ref)).

**Default `steady_state_solver` is `:auto`:** [`AbstractTraceData`](@ref) uses **`:augmented`** (AD-friendly steady state for [`loglikelihood_ad`](@ref)); other data use **`:default`**.

**Default `nuts_gradient` is `:auto`:** [`AbstractTraceData`](@ref) uses **`:finite`** (central differences; Zygote through the full trace HMM often overflows the compiler stack). RNA / count data use **`:Zygote`**. Pass **`nuts_gradient=:ForwardDiff`** to use forward-mode AD ([`NUTSOptions`](@ref); often good when only **few** parameters are fit). Pass `nuts_gradient=:Zygote` to try reverse-mode AD on traces.

**Traces and NUTS:** For live traces—including **coupled** `tracejoint` models—**`gradient=:finite`** is the **supported** NUTS choice: each log-density evaluation uses central differences on [`loglikelihood_ad`](@ref), so cost scales with the number of fitted parameters, not with “embarrassment.” Reverse-mode Zygote on the full coupled HMM (see [`benchmark_trace_zygote_gradient`](@ref)) is a separate, experimental micro-benchmark and can allocate huge memory or hit compiler `_pullback_generator` limits; it is **not** required for NUTS, and NUTS does **not** default to it for traces.

# Keywords
- `nchains`, `n_samples`, `n_adapts`, `nuts_δ`, `nuts_gradient`, `nuts_fd_ε`, `steady_state_solver` (`:auto` / `:default` / `:augmented`), `nuts_seed`, `ad_likelihood`: NUTS / likelihood options.
- `verbose`, `progress`: passed to [`NUTSOptions`](@ref) → [`AdvancedHMC.sample`](@ref). With **`nchains > 1`**, `progress` is forced to **`false`** on workers ([`Distributed.@spawn`](@ref) cannot reliably run AdvancedHMC’s progress machinery across processes). Use **`report_chains`** for parallel feedback, or run a **single chain** (`nchains=1`) if you want `progress=true`.
- `report_chains` (default `true`): log on the **main** process after each `fetch` (chain index order; use when runs are long).

# Returns
`NamedTuple` including `elapsed_sec`, merged `fits`, `stats`, `measures`, `chain_fits`, `nuts_results`.
"""
function benchmark_inference_run_nuts_parallel(
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel;
    nchains::Int=8,
    n_samples::Int=150,
    n_adapts::Int=80,
    nuts_δ::Float64=0.8,
    nuts_gradient::Symbol=:auto,
    nuts_fd_ε::Float64=1e-4,
    steady_state_solver::Symbol=:auto,
    nuts_seed::Int=2025,
    ad_likelihood::Union{Nothing,Bool}=nothing,
    verbose::Bool=false,
    progress::Bool=false,
    report_chains::Bool=true,
)
    ss = if steady_state_solver === :auto
        :augmented  # Use augmented solver for AD compatibility with NUTS
    else
        steady_state_solver
    end
    ng = if nuts_gradient === :auto
        data isa AbstractTraceData ? :finite : :Zygote
    else
        nuts_gradient
    end
    ng in (:Zygote, :finite, :ForwardDiff) ||
        throw(ArgumentError("nuts_gradient must be :auto, :Zygote, :finite, or :ForwardDiff, got $(repr(nuts_gradient))"))
    nchains > 1 && benchmark_inference_ensure_workers(nchains)
    progress_eff = progress
    if nchains > 1 && progress
        @warn "benchmark_inference_run_nuts_parallel: progress=true is not supported with nchains > 1 (Distributed); using progress=false on workers. Use report_chains=true or nchains=1."
        progress_eff = false
    end
    opts = NUTSOptions(;
        n_samples=n_samples,
        n_adapts=n_adapts,
        δ=nuts_δ,
        gradient=ng,
        fd_ε=nuts_fd_ε,
        verbose=verbose,
        progress=progress_eff,
    )
    rngs = [MersenneTwister(nuts_seed + k) for k in 1:nchains]
    report_chains && @info "benchmark_inference_run_nuts_parallel: spawning $(nchains) NUTS chains (n_samples=$(n_samples), n_adapts=$(n_adapts), steady_state_solver=$(ss), gradient=$(ng), nuts_fd_ε=$(nuts_fd_ε) (used if gradient is :finite))..."
    elapsed_sec = @elapsed begin
        results = if nchains == 1
            # Local run: avoids Distributed @spawn (ForwardDiff / progress / errors deserialize badly on workers).
            r = run_nuts_fit(data, model, opts; rng=rngs[1], steady_state_solver=ss, ad_likelihood=ad_likelihood)
            report_chains && @info "benchmark_inference_run_nuts_parallel: chain 1/1 finished"
            Any[r]
        else
            futures = map(1:nchains) do k
                @spawn run_nuts_fit(data, model, opts; rng=rngs[k], steady_state_solver=ss, ad_likelihood=ad_likelihood)
            end
            res = Vector{Any}(undef, nchains)
            for k in 1:nchains
                v = fetch(futures[k])
                v isa Distributed.RemoteException && throw(v)
                res[k] = v
                report_chains && @info "benchmark_inference_run_nuts_parallel: chain $(k)/$(nchains) finished"
            end
            res
        end
        for r in results
            r isa Distributed.RemoteException && throw(r)
            r isa Exception && throw(r)
            r isa Tuple || throw(ErrorException("benchmark_inference_run_nuts_parallel: expected each chain result to be a tuple from run_nuts_fit, got $(typeof(r))"))
        end
        fits_vec = [r[1] for r in results]
        merged = merge_fit(fits_vec)
        stats = compute_stats(merged.param, model)
        rhat = vec(compute_rhat(fits_vec))
        ess_pooled, geweke, mcse = compute_measures(fits_vec)
        waic = compute_waic(merged.lppd, merged.pwaic, data)
        measures = Measures(waic, rhat, ess_pooled, geweke, mcse)
    end
    ess_nt = _benchmark_ess_summary(measures, elapsed_sec)
    return (; elapsed_sec, fits=merged, stats, measures, chain_fits=fits_vec, nuts_results=results, ess_nt...)
end

function benchmark_inference_run_nuts_parallel(scenario::NamedTuple; kwargs...)
    hasproperty(scenario, :data) && hasproperty(scenario, :model) ||
        throw(ArgumentError("scenario must be a NamedTuple with :data and :model"))
    benchmark_inference_run_nuts_parallel(scenario.data, scenario.model; kwargs...)
end

"""
Longest per-trial frame count (for Zygote preflight).

Supports:
- `AbstractVector` of per-trial arrays (each trial `AbstractMatrix` or `AbstractVector`);
- `Tuple` with first component that `Vector` of per-trial arrays (**coupled `tracejoint`**: `(trials, background, weights, nframes)`);
- `Tuple` whose first component is a single `AbstractMatrix` (one stacked trial).
"""
_benchmark_trace_longest_timesteps(data::AbstractTraceData) = longest_trace_timesteps(data)

function _benchmark_zygote_nested_messages(err)::String
    parts = String[]
    e = err
    for _ in 1:6
        push!(parts, sprint(showerror, e))
        e isa TaskFailedException || break
        te = e.task.exception
        te === nothing && break
        e = te
    end
    return join(parts, "\n")
end

"""True if `err` looks like Zygote / compiler failure on a huge differentiated graph (not necessarily `StackOverflowError`)."""
function _benchmark_zygote_compiler_failure(err)::Bool
    err isa StackOverflowError && return true
    msg = _benchmark_zygote_nested_messages(err)
    lmsg = lowercase(msg)
    occursin("stack overflow", lmsg) && return true
    occursin("_pullback_generator", msg) && return true
    occursin("internal error", lmsg) && return true
    occursin("type inference", lmsg) && occursin("stack overflow", lmsg) && return true
    occursin("ntuple{", lmsg) && occursin("datatype", lmsg) && occursin("stack overflow", lmsg) && return true
    return false
end

"""
    benchmark_trace_zygote_gradient(data, model; kwargs...)
    benchmark_trace_zygote_gradient(scenario; kwargs...)

Time **Zygote** reverse-mode gradients of [`loglikelihood_ad`](@ref) at `get_param(model)` for trace
[`AbstractTraceData`](@ref), comparing **no** HMM checkpointing (`hmm_checkpoint_steps=nothing`) vs a
given chunk size. Uses the keyword `hmm_checkpoint_steps` on [`loglikelihood_ad`](@ref) (no global
[`set_hmm_zygote_checkpoint_steps!`](@ref)).

!!! warning "Long traces and Zygote"
    Reverse-mode AD through the **full** unrolled trace HMM often triggers **compiler stack overflow**
    (`InternalError` / `_pullback_generator` / huge `NTuple{…}`) when frame counts are in the thousands
    (Julia may print `Internal error` lines to stderr even when the run eventually finishes).
    By default this function **refuses** to run when the inferred longest trial exceeds `max_timesteps` (see below).
    [`benchmark_inference_trace_gr2r2`](@ref) defaults to **`totaltime=4000`** — that is usually **too long** for this
    benchmark; pass e.g. `totaltime=500.0, ntrials=2`, or use [`benchmark_trace_forwarddiff_gradient`](@ref) /
    [`benchmark_trace_finitediff_gradient`](@ref), or `skip_no_checkpoint=true`.

!!! danger "Coupled `tracejoint` models"
    If [`hastrait`](@ref)`(model, :coupling)` (e.g. [`benchmark_inference_trace_coupled_3x3`](@ref)), this helper **refuses**
    to run unless **`allow_coupled_zygote=true`**. Reverse-mode AD on the coupled full HMM can allocate **hundreds of GB**
    and is unsuitable for routine benchmarks. For **NUTS** on traces, [`benchmark_inference_run_nuts_parallel`](@ref) defaults
    to **`gradient=:finite`** (via `nuts_gradient=:auto`); that path does **not** use this Zygote helper.

Intended for **developer** micro-benchmarks (wall time, memory trade-offs on modest traces), not CI.
For timing gradients of a **small subset** of parameters (e.g. three), use
[`benchmark_trace_forwarddiff_gradient`](@ref) and [`benchmark_trace_finitediff_gradient`](@ref) (both support coupled
traces; large coupled state spaces can make ForwardDiff slow).

# Keywords
- `steady_state_solver` (default `:augmented`): passed to [`loglikelihood_ad`](@ref).
- `checkpoint_chunks` (default `32`): positive `Int` for the timed “with checkpoint” run; set to `nothing`
  to only report `time_no_checkpoint`.
- `nruns` (default `3`): each timing is the **minimum** of this many `@elapsed` trials (after one warmup each).
- `warmup` (default `true`): compile with one extra gradient per mode before timing.
- `max_timesteps` (default `2048`): if not `nothing`, require longest per-trial trace length (from `data.trace`)
  ≤ this, or throw `ArgumentError`. Pass `max_timesteps=nothing` to skip the check (may still hit compiler limits).
- `skip_no_checkpoint` (default `false`): if `true`, do **not** time the uncheckpointed gradient (the usual source
  of compiler `_pullback_generator` / `NTuple{…}` overflows on long traces); only time `hmm_checkpoint_steps=checkpoint_chunks`.
  Requires `checkpoint_chunks > 0`.
- `catch_zygote_failure` (default `true`): on Zygote/compiler stack overflow, return `success=false` and an
  `error` string instead of rethrowing. Set to `false` to debug with a full stack trace.
- `allow_coupled_zygote` (default `false`): if `hastrait(model, :coupling)`, throw unless `true` (see warning above).

# Returns
`NamedTuple` with `success::Bool`, `time_no_checkpoint`, `time_with_checkpoint` (or `nothing`), `ratio`,
`θlen`, `longest_timesteps` (or `nothing`), `steady_state_solver`, `checkpoint_chunks`. On failure,
`time_*` and `ratio` are `nothing` and `error` contains a short message.

# Example
```julia
using StochasticGene
scen = benchmark_inference_trace_gr2r2(totaltime=500.0, ntrials=2, seed=1)
benchmark_trace_zygote_gradient(scen; checkpoint_chunks=32, nruns=2)
```
"""
function benchmark_trace_zygote_gradient(
    data::AbstractTraceData,
    model::AbstractGRSMmodel;
    steady_state_solver::Symbol=:augmented,
    checkpoint_chunks::Union{Nothing,Integer}=32,
    nruns::Int=3,
    warmup::Bool=true,
    max_timesteps::Union{Nothing,Int}=2048,
    skip_no_checkpoint::Bool=false,
    catch_zygote_failure::Bool=true,
    allow_coupled_zygote::Bool=false,
)
    nruns < 1 && throw(ArgumentError("nruns must be >= 1"))
    max_timesteps !== nothing && max_timesteps < 0 &&
        throw(ArgumentError("max_timesteps must be non-negative or nothing"))
    skip_no_checkpoint && (checkpoint_chunks === nothing || checkpoint_chunks <= 0) &&
        throw(ArgumentError("skip_no_checkpoint=true requires checkpoint_chunks > 0"))
    if hastrait(model, :coupling) && !allow_coupled_zygote
        throw(ArgumentError(
            "benchmark_trace_zygote_gradient: coupled joint-trace model (hastrait(model, :coupling)); " *
            "Zygote reverse-mode on this graph can allocate hundreds of GB and is disabled by default. " *
            "Use benchmark_trace_forwarddiff_gradient / benchmark_trace_finitediff_gradient for subset gradients; " *
            "or pass allow_coupled_zygote=true only if you accept the memory and time cost.",
        ))
    end
    θ = Vector{Float64}(get_param(model))
    L = _benchmark_trace_longest_timesteps(data)
    if max_timesteps !== nothing && L !== nothing && L > max_timesteps
        throw(ArgumentError(
            "benchmark_trace_zygote_gradient: longest per-trial trace has $L frames (> max_timesteps=$max_timesteps). " *
            "Zygote reverse-mode on long trace HMMs usually overflows the compiler; use a shorter scenario " *
            "(e.g. smaller totaltime/ntrials), pass max_timesteps=nothing to skip this check at your own risk, " *
            "use skip_no_checkpoint=true with checkpoint_chunks>0 to time only the checkpointed path, " *
            "or use benchmark_trace_forwarddiff_gradient / benchmark_trace_finitediff_gradient.",
        ))
    end

    function mintime(hmm_checkpoint_steps)
        f = η -> loglikelihood_ad(
            η, data, model;
            steady_state_solver=steady_state_solver,
            hmm_stack=HMM_STACK_AD,
            hmm_checkpoint_steps=hmm_checkpoint_steps,
        )[1]
        warmup && Zygote.gradient(f, θ)
        t = Inf
        for _ in 1:nruns
            ti = @elapsed Zygote.gradient(f, θ)
            t = min(t, ti)
        end
        t
    end

    if skip_no_checkpoint
        local t1s
        try
            t1s = mintime(Int(checkpoint_chunks))
        catch err
            if catch_zygote_failure && _benchmark_zygote_compiler_failure(err)
                return (;
                    success=false,
                    time_no_checkpoint=nothing,
                    time_with_checkpoint=nothing,
                    ratio=nothing,
                    θlen=length(θ),
                    longest_timesteps=L,
                    steady_state_solver,
                    checkpoint_chunks=Int(checkpoint_chunks),
                    error=_benchmark_zygote_nested_messages(err),
                )
            end
            rethrow()
        end
        return (;
            success=true,
            time_no_checkpoint=nothing,
            time_with_checkpoint=t1s,
            ratio=nothing,
            θlen=length(θ),
            longest_timesteps=L,
            steady_state_solver,
            checkpoint_chunks=Int(checkpoint_chunks),
        )
    end

    local t0
    try
        t0 = mintime(nothing)
    catch err
        if catch_zygote_failure && _benchmark_zygote_compiler_failure(err)
            return (;
                success=false,
                time_no_checkpoint=nothing,
                time_with_checkpoint=nothing,
                ratio=nothing,
                θlen=length(θ),
                longest_timesteps=L,
                steady_state_solver,
                checkpoint_chunks,
                error=_benchmark_zygote_nested_messages(err),
            )
        end
        rethrow()
    end

    if checkpoint_chunks === nothing || checkpoint_chunks <= 0
        return (;
            success=true,
            time_no_checkpoint=t0,
            time_with_checkpoint=nothing,
            ratio=nothing,
            θlen=length(θ),
            longest_timesteps=L,
            steady_state_solver,
            checkpoint_chunks,
        )
    end

    local t1
    try
        t1 = mintime(Int(checkpoint_chunks))
    catch err
        if catch_zygote_failure && _benchmark_zygote_compiler_failure(err)
            return (;
                success=false,
                time_no_checkpoint=t0,
                time_with_checkpoint=nothing,
                ratio=nothing,
                θlen=length(θ),
                longest_timesteps=L,
                steady_state_solver,
                checkpoint_chunks=Int(checkpoint_chunks),
                error=_benchmark_zygote_nested_messages(err),
            )
        end
        rethrow()
    end
    return (;
        success=true,
        time_no_checkpoint=t0,
        time_with_checkpoint=t1,
        ratio=t1 / t0,
        θlen=length(θ),
        longest_timesteps=L,
        steady_state_solver,
        checkpoint_chunks=Int(checkpoint_chunks),
    )
end

function benchmark_trace_zygote_gradient(scenario::NamedTuple; kwargs...)
    hasproperty(scenario, :data) && hasproperty(scenario, :model) ||
        throw(ArgumentError("scenario must be a NamedTuple with :data and :model"))
    benchmark_trace_zygote_gradient(scenario.data, scenario.model; kwargs...)
end

function _benchmark_trace_resolve_param_indices(θlen::Int, param_indices::Union{Nothing,AbstractVector{<:Integer}})
    if param_indices === nothing
        k = min(3, θlen)
        k < 1 && throw(ArgumentError("parameter vector is empty"))
        return collect(1:k)
    end
    idx = collect(Int, param_indices)
    isempty(idx) && throw(ArgumentError("param_indices must be non-empty"))
    all(i -> 1 <= i <= θlen, idx) || throw(ArgumentError("param_indices must lie in 1:$θlen"))
    return idx
end

"""One central-difference gradient vector for `θ` at `idx` (same convention as [`benchmark_trace_finitediff_gradient`](@ref))."""
function _benchmark_trace_central_grad_vector(
    θ::Vector{Float64},
    idx::Vector{Int},
    data::AbstractTraceData,
    model::AbstractGRSMmodel,
    fd_ε::Float64,
    steady_state_solver::Symbol,
    hmm_checkpoint_steps::Union{Nothing,Integer},
)
    k = length(idx)
    g = Vector{Float64}(undef, k)
    for j in 1:k
        i = idx[j]
        θp = copy(θ)
        θm = copy(θ)
        θp[i] += fd_ε
        θm[i] -= fd_ε
        lp = loglikelihood_ad(
            θp, data, model;
            steady_state_solver=steady_state_solver,
            hmm_stack=HMM_STACK_AD,
            hmm_checkpoint_steps=hmm_checkpoint_steps,
        )[1]
        lm = loglikelihood_ad(
            θm, data, model;
            steady_state_solver=steady_state_solver,
            hmm_stack=HMM_STACK_AD,
            hmm_checkpoint_steps=hmm_checkpoint_steps,
        )[1]
        g[j] = (lp - lm) / (2fd_ε)
    end
    return g
end

function _benchmark_trace_ll_subparams(
    u::AbstractVector,
    θ::Vector{Float64},
    idx::Vector{Int},
    data::AbstractTraceData,
    model::AbstractGRSMmodel;
    steady_state_solver::Symbol,
    hmm_checkpoint_steps::Union{Nothing,Integer},
)
    seed = first(u)
    θbase = [θi + zero(seed) for θi in θ]
    θc = [
        begin
            j = findfirst(==(i), idx)
            j === nothing ? θbase[i] : u[j]
        end
        for i in eachindex(θbase)
    ]
    return loglikelihood_ad(
        θc, data, model;
        steady_state_solver=steady_state_solver,
        hmm_stack=HMM_STACK_AD,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
    )[1]
end

function _benchmark_trace_fwd_vs_fd_agreement(
    θ::Vector{Float64},
    idx::Vector{Int},
    data::AbstractTraceData,
    model::AbstractGRSMmodel,
    fd_ε::Float64,
    steady_state_solver::Symbol,
    hmm_checkpoint_steps::Union{Nothing,Integer},
)
    u0 = θ[idx]
    function ll_u(u)
        _benchmark_trace_ll_subparams(u, θ, idx, data, model; steady_state_solver=steady_state_solver, hmm_checkpoint_steps=hmm_checkpoint_steps)
    end
    g_fwd = ForwardDiff.gradient(ll_u, u0)
    g_fd = _benchmark_trace_central_grad_vector(
        θ, idx, data, model, fd_ε, steady_state_solver, hmm_checkpoint_steps,
    )
    mad = maximum(abs, g_fwd .- g_fd)
    denom = max(maximum(abs, g_fwd), maximum(abs, g_fd), 1e-12)
    return (; gradient_max_abs_diff=mad, gradient_max_rel_diff=mad / denom)
end

function _benchmark_trace_zygote_vs_fd_agreement(
    θ::Vector{Float64},
    idx::Vector{Int},
    data::AbstractTraceData,
    model::AbstractGRSMmodel,
    fd_ε::Float64,
    steady_state_solver::Symbol,
    hmm_checkpoint_steps::Union{Nothing,Integer},
)
    u0 = θ[idx]
    function ll_u(u)
        _benchmark_trace_ll_subparams(u, θ, idx, data, model; steady_state_solver=steady_state_solver, hmm_checkpoint_steps=hmm_checkpoint_steps)
    end
    gz = Zygote.gradient(ll_u, u0)
    g_z = gz === nothing ? nothing : first(gz)
    g_z === nothing && error("Zygote.gradient returned nothing for subset log-likelihood")
    g_fd = _benchmark_trace_central_grad_vector(
        θ, idx, data, model, fd_ε, steady_state_solver, hmm_checkpoint_steps,
    )
    mad = maximum(abs, g_z .- g_fd)
    denom = max(maximum(abs, g_z), maximum(abs, g_fd), 1e-12)
    return (; gradient_max_abs_diff=mad, gradient_max_rel_diff=mad / denom)
end

"""
    benchmark_trace_forwarddiff_gradient(data, model; kwargs...)
    benchmark_trace_forwarddiff_gradient(scenario; kwargs...)

Wall-clock time for a **full forward-mode gradient** of [`loglikelihood_ad`](@ref) with respect to a
**subset** of transformed parameters (default: the first **three** entries of `get_param(model)`), holding
the rest fixed. This mirrors the use case for NUTS with [`NUTSOptions`](@ref) `gradient=:ForwardDiff` when
only a few coordinates matter or dimension is small.

Uses [`ForwardDiff.gradient`](https://juliadiff.org/ForwardDiff.jl/stable/user/advanced/#ForwardDiff.gradient)
on `u ↦ loglikelihood_ad(θ with u embedded at param_indices, …)`.

# Keywords
- `param_indices`: indices into `θ` (1-based); default `nothing` → `1:min(3, length(θ))`.
- `steady_state_solver` (default `:augmented`), `hmm_checkpoint_steps` (default `nothing`): passed to [`loglikelihood_ad`](@ref).
- `fd_ε`: accepted for compatibility with shared benchmark keyword bundles; **ignored** by the ForwardDiff path.
- `nruns`, `warmup`: same meaning as [`benchmark_trace_zygote_gradient`](@ref).

# Returns
`NamedTuple` with `method` (`:forwarddiff`), `time_sec`, **`gradient`** (`Vector{Float64}`, last timed evaluation),
`param_indices`, `nparams`, `θlen`, `steady_state_solver`, `hmm_checkpoint_steps`.

# Example
```julia
using StochasticGene
scen = benchmark_inference_trace_gr2r2(totaltime=400.0, ntrials=4, seed=1)
benchmark_trace_forwarddiff_gradient(scen; param_indices=[1, 2, 4], nruns=2)
```
"""
function benchmark_trace_forwarddiff_gradient(
    data::AbstractTraceData,
    model::AbstractGRSMmodel;
    param_indices::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    steady_state_solver::Symbol=:augmented,
    hmm_checkpoint_steps::Union{Nothing,Integer}=nothing,
    fd_ε::Union{Nothing,Float64}=nothing,
    nruns::Int=3,
    warmup::Bool=true,
)
    nruns < 1 && throw(ArgumentError("nruns must be >= 1"))
    fd_ε !== nothing && fd_ε <= 0 && throw(ArgumentError("fd_ε must be positive when provided"))
    θ = Vector{Float64}(get_param(model))
    θlen = length(θ)
    idx = _benchmark_trace_resolve_param_indices(θlen, param_indices)
    u0 = θ[idx]
    function ll_u(u)
        _benchmark_trace_ll_subparams(u, θ, idx, data, model; steady_state_solver=steady_state_solver, hmm_checkpoint_steps=hmm_checkpoint_steps)
    end
    warmup && ForwardDiff.gradient(ll_u, u0)
    t = Inf
    local g_out
    for _ in 1:nruns
        ti = @elapsed g_out = ForwardDiff.gradient(ll_u, u0)
        t = min(t, ti)
    end
    return (;
        method=:forwarddiff,
        time_sec=t,
        gradient=collect(Float64, g_out),
        param_indices=idx,
        nparams=length(idx),
        θlen,
        steady_state_solver,
        hmm_checkpoint_steps,
    )
end

function benchmark_trace_forwarddiff_gradient(scenario::NamedTuple; kwargs...)
    hasproperty(scenario, :data) && hasproperty(scenario, :model) ||
        throw(ArgumentError("scenario must be a NamedTuple with :data and :model"))
    benchmark_trace_forwarddiff_gradient(scenario.data, scenario.model; kwargs...)
end

"""
    benchmark_trace_finitediff_gradient(data, model; kwargs...)
    benchmark_trace_finitediff_gradient(scenario; kwargs...)

Wall-clock time for a **central finite-difference gradient** of [`loglikelihood_ad`](@ref) with respect to a
subset of parameters (default: first **three**), matching NUTS `gradient=:finite` (central difference per
coordinate). Each gradient uses **``2k``** likelihood evaluations for `k` parameters.

# Keywords
- `param_indices`: default `nothing` → `1:min(3, length(θ))`.
- `fd_ε` (default `1e-4`): same order of magnitude as [`NUTSOptions`](@ref) `fd_ε`.
- `steady_state_solver`, `hmm_checkpoint_steps`: passed to [`loglikelihood_ad`](@ref).
- `nruns`, `warmup`: same as [`benchmark_trace_forwarddiff_gradient`](@ref).

# Returns
`NamedTuple` with `method` (`:finitediff`), `time_sec`, **`gradient`** (`Vector{Float64}`, last timed evaluation),
`param_indices`, `nparams`, `n_ll_evals_per_grad` (= `2 * nparams`), `fd_ε`, `θlen`, etc.

# Example
```julia
using StochasticGene
scen = benchmark_inference_trace_gr2r2(totaltime=400.0, ntrials=4, seed=1)
benchmark_trace_finitediff_gradient(scen; param_indices=collect(1:3), fd_ε=1e-4, nruns=2)
```
"""
function benchmark_trace_finitediff_gradient(
    data::AbstractTraceData,
    model::AbstractGRSMmodel;
    param_indices::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    steady_state_solver::Symbol=:augmented,
    hmm_checkpoint_steps::Union{Nothing,Integer}=nothing,
    fd_ε::Float64=1e-4,
    nruns::Int=3,
    warmup::Bool=true,
)
    nruns < 1 && throw(ArgumentError("nruns must be >= 1"))
    fd_ε > 0 || throw(ArgumentError("fd_ε must be positive"))
    θ = Vector{Float64}(get_param(model))
    θlen = length(θ)
    idx = _benchmark_trace_resolve_param_indices(θlen, param_indices)
    k = length(idx)
    function fd_grad!(g::Vector{Float64}, θ0::Vector{Float64})
        @assert length(g) == k
        g .= _benchmark_trace_central_grad_vector(
            θ0, idx, data, model, fd_ε, steady_state_solver, hmm_checkpoint_steps,
        )
        return g
    end
    g = Vector{Float64}(undef, k)
    warmup && fd_grad!(g, θ)
    t = Inf
    for _ in 1:nruns
        ti = @elapsed fd_grad!(g, θ)
        t = min(t, ti)
    end
    return (;
        method=:finitediff,
        time_sec=t,
        gradient=copy(g),
        param_indices=idx,
        nparams=k,
        n_ll_evals_per_grad=2k,
        fd_ε,
        θlen,
        steady_state_solver,
        hmm_checkpoint_steps,
    )
end

function benchmark_trace_finitediff_gradient(scenario::NamedTuple; kwargs...)
    hasproperty(scenario, :data) && hasproperty(scenario, :model) ||
        throw(ArgumentError("scenario must be a NamedTuple with :data and :model"))
    benchmark_trace_finitediff_gradient(scenario.data, scenario.model; kwargs...)
end

"""
    benchmark_trace_zygote_subset_gradient(data, model; kwargs...)
    benchmark_trace_zygote_subset_gradient(scenario; kwargs...)

Wall-clock time for **Zygote** reverse-mode gradient of [`loglikelihood_ad`](@ref) on the same **subset** of
transformed parameters as [`benchmark_trace_forwarddiff_gradient`](@ref) / [`benchmark_trace_finitediff_gradient`](@ref)
(`u ↦` likelihood with `u` embedded at `param_indices`).

Uses [`Zygote.gradient`](https://fluxml.ai/Zygote.jl/stable/) on that scalar map. Coupled joint traces can be **slow**
or memory-heavy; use [`hmm_checkpoint_steps`](@ref) and per-trace checkpointing in the HMM stack.

# Keywords
Same as [`benchmark_trace_forwarddiff_gradient`](@ref) (`param_indices`, `steady_state_solver`, `hmm_checkpoint_steps`, `nruns`, `warmup`).
`fd_ε` is also accepted for compatibility with shared keyword bundles and ignored on the Zygote path.

# Returns
`NamedTuple` with `method` (`:zygote`), `time_sec`, **`gradient`**, `param_indices`, `nparams`, `θlen`, `steady_state_solver`, `hmm_checkpoint_steps`.
"""
function benchmark_trace_zygote_subset_gradient(
    data::AbstractTraceData,
    model::AbstractGRSMmodel;
    param_indices::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    steady_state_solver::Symbol=:augmented,
    hmm_checkpoint_steps::Union{Nothing,Integer}=nothing,
    fd_ε::Union{Nothing,Float64}=nothing,
    nruns::Int=3,
    warmup::Bool=true,
)
    nruns < 1 && throw(ArgumentError("nruns must be >= 1"))
    fd_ε !== nothing && fd_ε <= 0 && throw(ArgumentError("fd_ε must be positive when provided"))
    θ = Vector{Float64}(get_param(model))
    θlen = length(θ)
    idx = _benchmark_trace_resolve_param_indices(θlen, param_indices)
    u0 = θ[idx]
    function ll_u(u)
        _benchmark_trace_ll_subparams(u, θ, idx, data, model; steady_state_solver=steady_state_solver, hmm_checkpoint_steps=hmm_checkpoint_steps)
    end
    warmup && Zygote.gradient(ll_u, u0)
    t = Inf
    local gz
    for _ in 1:nruns
        ti = @elapsed gz = Zygote.gradient(ll_u, u0)
        t = min(t, ti)
    end
    g1 = gz === nothing ? nothing : first(gz)
    g1 === nothing && throw(ErrorException("Zygote.gradient returned nothing for trace subset likelihood"))
    return (;
        method=:zygote,
        time_sec=t,
        gradient=collect(Float64, g1),
        param_indices=idx,
        nparams=length(idx),
        θlen,
        steady_state_solver,
        hmm_checkpoint_steps,
    )
end

function benchmark_trace_zygote_subset_gradient(scenario::NamedTuple; kwargs...)
    hasproperty(scenario, :data) && hasproperty(scenario, :model) ||
        throw(ArgumentError("scenario must be a NamedTuple with :data and :model"))
    benchmark_trace_zygote_subset_gradient(scenario.data, scenario.model; kwargs...)
end

function _benchmark_trace_pairwise_grad_diff(a::AbstractVector{Float64}, b::AbstractVector{Float64})
    d = a .- b
    mad = maximum(abs, d)
    denom = max(maximum(abs, a), maximum(abs, b), 1e-12)
    return (; max_abs_diff=mad, max_rel_diff=mad / denom)
end

"""
    compare_trace_subset_gradient_benchmarks(r_forwarddiff, r_finitediff, r_zygote)

Compare outputs from [`benchmark_trace_forwarddiff_gradient`](@ref), [`benchmark_trace_finitediff_gradient`](@ref),
and [`benchmark_trace_zygote_subset_gradient`](@ref) run with the **same** `data`, `model`, and keyword options.

# Arguments
- `r_forwarddiff`, `r_finitediff`, `r_zygote`: `NamedTuple`s from the three benchmark functions (must share `param_indices`, `nparams`, `θlen`).

# Returns
`NamedTuple` with `compare` (`:trace_subset_gradient_three_way`), wall-clock times, `fastest_method` (`:forwarddiff`, `:finitediff`, or `:zygote`),
pairwise `max_abs_diff` / `max_rel_diff` between the three gradient vectors, and shared metadata (`param_indices`, `fd_ε`, etc.).

# Example
```julia
using StochasticGene
scen = benchmark_scenario_coupled_3x3_10traces_220frames(seed=1)
rf = benchmark_trace_forwarddiff_gradient(scen; param_indices=[1], nruns=2, hmm_checkpoint_steps=32)
rfd = benchmark_trace_finitediff_gradient(scen; param_indices=[1], nruns=2, hmm_checkpoint_steps=32, fd_ε=1e-4)
rz = benchmark_trace_zygote_subset_gradient(scen; param_indices=[1], nruns=2, hmm_checkpoint_steps=32)
compare_trace_subset_gradient_benchmarks(rf, rfd, rz)
```
"""
function compare_trace_subset_gradient_benchmarks(
    r_forwarddiff::NamedTuple,
    r_finitediff::NamedTuple,
    r_zygote::NamedTuple,
)
    for k in (:param_indices, :nparams, :θlen)
        a = getfield(r_forwarddiff, k)
        getfield(r_finitediff, k) == a && getfield(r_zygote, k) == a ||
            throw(ArgumentError("compare_trace_subset_gradient_benchmarks: mismatch in $(k) across results"))
    end
    gf = r_forwarddiff.gradient
    gfd = r_finitediff.gradient
    gz = r_zygote.gradient
    length(gf) == length(gfd) == length(gz) || throw(ArgumentError("compare_trace_subset_gradient_benchmarks: gradient length mismatch"))
    d_ffd = _benchmark_trace_pairwise_grad_diff(gf, gfd)
    d_fz = _benchmark_trace_pairwise_grad_diff(gf, gz)
    d_fdz = _benchmark_trace_pairwise_grad_diff(gfd, gz)
    ts = [r_forwarddiff.time_sec, r_finitediff.time_sec, r_zygote.time_sec]
    fastest_method = [:forwarddiff, :finitediff, :zygote][argmin(ts)]
    return (;
        compare=:trace_subset_gradient_three_way,
        time_forwarddiff_sec=r_forwarddiff.time_sec,
        time_finitediff_sec=r_finitediff.time_sec,
        time_zygote_sec=r_zygote.time_sec,
        fastest_method,
        forwarddiff_vs_finitediff_max_abs_diff=d_ffd.max_abs_diff,
        forwarddiff_vs_finitediff_max_rel_diff=d_ffd.max_rel_diff,
        forwarddiff_vs_zygote_max_abs_diff=d_fz.max_abs_diff,
        forwarddiff_vs_zygote_max_rel_diff=d_fz.max_rel_diff,
        finitediff_vs_zygote_max_abs_diff=d_fdz.max_abs_diff,
        finitediff_vs_zygote_max_rel_diff=d_fdz.max_rel_diff,
        param_indices=r_forwarddiff.param_indices,
        nparams=r_forwarddiff.nparams,
        θlen=r_forwarddiff.θlen,
        fd_ε=get(r_finitediff, :fd_ε, nothing),
        steady_state_solver=r_forwarddiff.steady_state_solver,
        hmm_checkpoint_steps=r_forwarddiff.hmm_checkpoint_steps,
    )
end

"""
    benchmark_trace_compare_forwarddiff_vs_finitediff(data, model; kwargs...)
    benchmark_trace_compare_forwarddiff_vs_finitediff(scenario; kwargs...)

Compare wall-clock time for a **subset gradient** of [`loglikelihood_ad`](@ref) using **forward-mode**
([`benchmark_trace_forwarddiff_gradient`](@ref)) vs **central finite differences**
([`benchmark_trace_finitediff_gradient`](@ref)), with identical `param_indices`, `steady_state_solver`,
`hmm_checkpoint_steps`, `nruns`, and `warmup`. Use this to choose between NUTS `gradient=:ForwardDiff` and
`gradient=:finite` for the same trace scenario and number of parameters.

Coupled joint traces are supported; large coupled state spaces may make ForwardDiff gradients slow.

# Keywords
- `param_indices`, `steady_state_solver`, `hmm_checkpoint_steps`, `nruns`, `warmup`: shared by both benchmarks.
- `fd_ε`: passed only to the finite-difference run (default `1e-4`).
- `check_gradients` (default `true`): after timing, compare an independent AD gradient to central FD on the same
  `u ↦ loglikelihood_ad(θ with u at param_indices, …)`:
  - **ForwardDiff vs FD** (`_benchmark_trace_fwd_vs_fd_agreement`) when `check_gradients=true`.
  - Optionally set **`check_zygote_vs_fd=true`** on **coupled** models to **also** compare **Zygote** vs FD (expensive,
    can OOM on long traces — same caveats as [`benchmark_trace_zygote_gradient`](@ref)); default `false`.
- `check_zygote_vs_fd` (default `false`): coupled models only; run Zygote vs central FD in addition to ForwardDiff vs FD.

# Returns
`NamedTuple` with:
- `compare` (`:forwarddiff_vs_finitediff`): tag for [`benchmark_inference_print_summary`](@ref)
- `time_forwarddiff_sec`: seconds for one ForwardDiff gradient
- `time_finitediff_sec`: seconds for one central finite-difference gradient
- `faster`: `:forwarddiff` or `:finite`
- `ratio_forwarddiff_over_finitediff`: `time_forwarddiff_sec / time_finitediff_sec`; **< 1** means ForwardDiff is faster, **> 1** means finite differences are faster
- `gradient_check_reference`: `:forwarddiff_vs_fd`, `:zygote_vs_fd` (extra coupled check), `:none` (checks off), `:failed`, or `:unavailable`
- `gradient_max_abs_diff`, `gradient_max_rel_diff`: agreement for the reference pair above (`nothing` when no check ran).
  If components are very small, rely on `gradient_max_abs_diff`; FD truncation error scales with `fd_ε`.
- `forwarddiff_unavailable_reason`: `nothing` or a short string when ForwardDiff was skipped
- `param_indices`, `nparams`, `n_ll_evals_per_grad`, `fd_ε`, `θlen`, `steady_state_solver`, `hmm_checkpoint_steps`

# Example
```julia
using StochasticGene
scen = benchmark_inference_trace_gr2r2(totaltime=200.0, ntrials=2, seed=1)
benchmark_trace_compare_forwarddiff_vs_finitediff(scen; param_indices=[1, 2], nruns=3, fd_ε=1e-4)
```

Wrapping the **first** call in `@time` in a fresh session often shows most time in **compilation**; the returned
`time_*_sec` fields are still minima over `nruns` after internal warmup (`warmup=true` by default). Repeat the
call (or run after precompilation) if you need `@time` without the one-time compile cost.
"""
function benchmark_trace_compare_forwarddiff_vs_finitediff(
    data::AbstractTraceData,
    model::AbstractGRSMmodel;
    param_indices::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    steady_state_solver::Symbol=:augmented,
    hmm_checkpoint_steps::Union{Nothing,Integer}=nothing,
    fd_ε::Float64=1e-4,
    nruns::Int=3,
    warmup::Bool=true,
    check_gradients::Bool=true,
    check_zygote_vs_fd::Bool=false,
)
    fd = benchmark_trace_finitediff_gradient(
        data, model;
        param_indices=param_indices,
        steady_state_solver=steady_state_solver,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
        fd_ε=fd_ε,
        nruns=nruns,
        warmup=warmup,
    )
    fw = benchmark_trace_forwarddiff_gradient(
        data, model;
        param_indices=param_indices,
        steady_state_solver=steady_state_solver,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
        nruns=nruns,
        warmup=warmup,
    )
    tf = fw.time_sec
    tfd = fd.time_sec
    agr = if check_gradients
        θ = Vector{Float64}(get_param(model))
        _benchmark_trace_fwd_vs_fd_agreement(
            θ, fd.param_indices, data, model, fd_ε, steady_state_solver, hmm_checkpoint_steps,
        )
    else
        (; gradient_max_abs_diff=nothing, gradient_max_rel_diff=nothing)
    end
    if check_gradients && hastrait(model, :coupling) && check_zygote_vs_fd
        θz = Vector{Float64}(get_param(model))
        try
            zagr = _benchmark_trace_zygote_vs_fd_agreement(
                θz, fd.param_indices, data, model, fd_ε, steady_state_solver, hmm_checkpoint_steps,
            )
            @info "benchmark_trace_compare (coupled): optional Zygote vs FD" gradient_max_abs_diff=zagr.gradient_max_abs_diff gradient_max_rel_diff=zagr.gradient_max_rel_diff
        catch e
            @warn "Zygote vs FD (optional coupled check) failed; set check_zygote_vs_fd=false to skip" exception=(e, catch_backtrace())
        end
    end
    return (;
        compare=:forwarddiff_vs_finitediff,
        time_forwarddiff_sec=tf,
        time_finitediff_sec=tfd,
        faster=tf < tfd ? :forwarddiff : :finite,
        ratio_forwarddiff_over_finitediff=tf / tfd,
        gradient_check_reference=check_gradients ? :forwarddiff_vs_fd : :none,
        gradient_max_abs_diff=agr.gradient_max_abs_diff,
        gradient_max_rel_diff=agr.gradient_max_rel_diff,
        forwarddiff_unavailable_reason=nothing,
        param_indices=fd.param_indices,
        nparams=fd.nparams,
        n_ll_evals_per_grad=fd.n_ll_evals_per_grad,
        fd_ε=fd.fd_ε,
        θlen=fd.θlen,
        steady_state_solver=fd.steady_state_solver,
        hmm_checkpoint_steps=fd.hmm_checkpoint_steps,
    )
end

function benchmark_trace_compare_forwarddiff_vs_finitediff(scenario::NamedTuple; kwargs...)
    hasproperty(scenario, :data) && hasproperty(scenario, :model) ||
        throw(ArgumentError("scenario must be a NamedTuple with :data and :model"))
    benchmark_trace_compare_forwarddiff_vs_finitediff(scenario.data, scenario.model; kwargs...)
end

"""
    benchmark_inference_run_advi(data, model; kwargs...)
    benchmark_inference_run_advi(scenario; kwargs...)

Mean-field ADVI benchmark runner with the **same comparison-facing payload shape** as the MH and NUTS
benchmark helpers: it returns `fits`, `stats`, `measures`, plus ADVI-specific diagnostics.
This makes scenario-based head-to-head comparisons straightforward.

For benchmark use, the default ADVI settings are intentionally **interactive** rather than
maximum-effort: `gradient=:finite` (implemented via `ForwardDiff` on `neg_elbo`), moderate
`maxiter`, and a bounded `time_limit`. This makes head-to-head scenario comparisons practical.
If you explicitly request `gradient=:Zygote` and that path fails on a scenario, this helper
retries once with `:finite` so the benchmark still completes.

# Keywords
- `rng`, `steady_state_solver`, `ad_likelihood`, `zygote_trace`: passed to [`run_advi_fit`](@ref).
- `advi_options`: [`ADVIOptions`](@ref) for the ADVI optimization. The benchmark default uses
  `gradient=:finite` for robust comparisons.

# Returns
`NamedTuple` including `elapsed_sec`, `fits`, `stats`, `measures`, `advi_info`, `neg_elbo_min`,
`gradient_used`, `zygote_fallback_used`, and the same ESS summary fields used by the MH/NUTS
benchmark runners (typically `NaN`/formal-only for ADVI).
"""
function benchmark_inference_run_advi(
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel;
    rng::Random.AbstractRNG=Random.default_rng(),
    steady_state_solver::Symbol=:auto,
    ad_likelihood::Union{Nothing,Bool}=nothing,
    zygote_trace::Bool=false,
    posterior_samples::Int=1000,
    advi_options::ADVIOptions=ADVIOptions(; maxiter=200, n_mc=4, σ_floor=1e-4, verbose=false, gradient=:finite, time_limit=60.0),
)
    ss = if steady_state_solver === :auto
        data isa AbstractTraceData ? :augmented : :default
    else
        steady_state_solver
    end
    opts_used = advi_options
    zygote_fallback_used = false
    local fits, stats, measures, advi_info
    elapsed_sec = @elapsed begin
        try
            fits, stats, measures, advi_info = run_advi_fit(
                data, model, opts_used;
                rng=rng,
                steady_state_solver=ss,
                ad_likelihood=ad_likelihood,
                zygote_trace=zygote_trace,
                posterior_samples=posterior_samples,
            )
        catch e
            if opts_used.gradient === :Zygote
                @warn "benchmark_inference_run_advi: ADVI with gradient=:Zygote failed on this scenario; retrying with gradient=:finite for a robust benchmark comparison."
                opts_used = ADVIOptions(
                    opts_used.maxiter,
                    opts_used.n_mc,
                    opts_used.σ_floor,
                    opts_used.init_s_raw,
                    opts_used.verbose,
                    :finite,
                    opts_used.time_limit,
                )
                zygote_fallback_used = true
                fits, stats, measures, advi_info = run_advi_fit(
                    data, model, opts_used;
                    rng=rng,
                    steady_state_solver=ss,
                    ad_likelihood=ad_likelihood,
                    zygote_trace=zygote_trace,
                    posterior_samples=posterior_samples,
                )
            else
                rethrow(e)
            end
        end
    end
    ess = _benchmark_ess_summary(measures, elapsed_sec)
    neg_elbo_min = Optim.minimum(advi_info.optimization)
    optimization_converged = Optim.converged(advi_info.optimization)
    optimizer_iterations = Optim.iterations(advi_info.optimization)
    posterior_mean = vec(stats.meanparam)
    posterior_std = vec(stats.stdparam)
    posterior_median = vec(stats.medparam)
    posterior_q = stats.qparam
    return (; elapsed_sec, fits, stats, measures, advi_info, neg_elbo_min, μ=advi_info.μ, σ=advi_info.σ, optimization=advi_info.optimization, initial_θ=advi_info.initial_θ, gradient_used=opts_used.gradient, zygote_fallback_used, steady_state_solver_used=ss, optimization_converged, optimizer_iterations, posterior_samples, posterior_mean, posterior_std, posterior_median, posterior_q, ess...)
end

function benchmark_inference_run_advi(scenario::NamedTuple; kwargs...)
    hasproperty(scenario, :data) && hasproperty(scenario, :model) ||
        throw(ArgumentError("scenario must be a NamedTuple with :data and :model"))
    benchmark_inference_run_advi(scenario.data, scenario.model; kwargs...)
end

"""
    benchmark_inference_compare_mh_nuts_advi(data, model; kwargs...)
    benchmark_inference_compare_mh_nuts_advi(scenario; kwargs...)

Run **MH**, **NUTS**, and **ADVI** on the same benchmark target so simulated scenarios from
[`benchmark_inference_simrna_small`](@ref), [`benchmark_inference_trace_gr2r2`](@ref), or related
helpers can be compared head-to-head without rebuilding `data` / `model`.

This is intentionally thin: it reuses [`benchmark_inference_run_mh`](@ref),
[`benchmark_inference_run_nuts_parallel`](@ref), and [`benchmark_inference_run_advi`](@ref)
rather than duplicating inference code.

# Keywords
- `run_mh`, `run_nuts`, `run_advi`: enable/disable individual methods.
- `mh_kwargs`, `nuts_kwargs`, `advi_kwargs`: `NamedTuple`s splatted into the corresponding
  benchmark runner.

# Returns
A `NamedTuple` with keys `mh`, `nuts`, `advi`, and `compare=:mh_nuts_advi`.
Each field is either the benchmark result `NamedTuple` from the underlying runner or `nothing`
if that method was disabled.

# Example
```julia
scen = benchmark_inference_simrna_small(seed=1)
out = benchmark_inference_compare_mh_nuts_advi(
    scen;
    mh_kwargs=(; nchains=1, samplesteps=2_000, warmupsteps=1_000, maxtime=30.0),
    nuts_kwargs=(; nchains=1, n_samples=50, n_adapts=50, report_chains=false),
    advi_kwargs=(; advi_options=ADVIOptions(maxiter=300, gradient=:finite, time_limit=60.0)),
)
benchmark_inference_print_summary(out)
```
"""
function benchmark_inference_compare_mh_nuts_advi(
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel;
    run_mh::Bool=true,
    run_nuts::Bool=true,
    run_advi::Bool=true,
    mh_kwargs::NamedTuple=(; nchains=1, samplesteps=10_000, warmupsteps=5_000, maxtime=600.0),
    nuts_kwargs::NamedTuple=(; nchains=1, n_samples=150, n_adapts=80, report_chains=true, verbose=false, progress=false),
    advi_kwargs::NamedTuple=(;),
)
    mh = run_mh ? benchmark_inference_run_mh(data, model; mh_kwargs...) : nothing
    nuts = run_nuts ? benchmark_inference_run_nuts_parallel(data, model; nuts_kwargs...) : nothing
    advi = run_advi ? benchmark_inference_run_advi(data, model; advi_kwargs...) : nothing
    return (; compare=:mh_nuts_advi, mh, nuts, advi)
end

function benchmark_inference_compare_mh_nuts_advi(scenario::NamedTuple; kwargs...)
    hasproperty(scenario, :data) && hasproperty(scenario, :model) ||
        throw(ArgumentError("scenario must be a NamedTuple with :data and :model"))
    out = benchmark_inference_compare_mh_nuts_advi(scenario.data, scenario.model; kwargs...)
    return hasproperty(scenario, :meta) ? merge(out, (; meta=scenario.meta)) : out
end

"""
    benchmark_inference_print_summary(result::NamedTuple)

Print a short summary for benchmark `NamedTuple` returns:

- [`benchmark_inference_run_mh`](@ref) / [`benchmark_inference_run_nuts_parallel`](@ref): `elapsed_sec` and, when present, ESS / R-hat lines.
- [`benchmark_inference_run_advi`](@ref): `elapsed_sec`, `neg_elbo_min`, and now `fits` / `stats` / `measures` for direct comparison with MH and NUTS.
- [`benchmark_trace_zygote_gradient`](@ref): `success`, `time_no_checkpoint`, `time_with_checkpoint`, `ratio`, `θlen`, `longest_timesteps`, and `error` on failure.
- [`benchmark_trace_forwarddiff_gradient`](@ref) / [`benchmark_trace_finitediff_gradient`](@ref): `time_sec`, `param_indices`, `nparams`, `θlen`, and finite-diff extras when present.
- [`benchmark_trace_compare_forwarddiff_vs_finitediff`](@ref): `time_forwarddiff_sec`, `time_finitediff_sec`, `faster`, `ratio_forwarddiff_over_finitediff`.
- [`compare_trace_subset_gradient_benchmarks`](@ref): `time_*_sec`, `fastest_method`, pairwise gradient diffs.
"""
function benchmark_inference_print_summary(result::NamedTuple)
    if get(result, :compare, nothing) === :mh_nuts_advi
        println("benchmark_inference_compare_mh_nuts_advi:")
        if get(result, :mh, nothing) !== nothing
            println("  MH elapsed_sec: ", result.mh.elapsed_sec)
            haskey(result.mh, :min_ess) && println("  MH min_ess / mean_ess / max_rhat: ", result.mh.min_ess, " / ", result.mh.mean_ess, " / ", result.mh.max_rhat)
            haskey(result.mh, :ess_per_sec_min) && println("  MH ESS/s (min / mean): ", result.mh.ess_per_sec_min, " / ", result.mh.ess_per_sec_mean)
        end
        if get(result, :nuts, nothing) !== nothing
            println("  NUTS elapsed_sec: ", result.nuts.elapsed_sec)
            haskey(result.nuts, :min_ess) && println("  NUTS min_ess / mean_ess / max_rhat: ", result.nuts.min_ess, " / ", result.nuts.mean_ess, " / ", result.nuts.max_rhat)
            haskey(result.nuts, :ess_per_sec_min) && println("  NUTS ESS/s (min / mean): ", result.nuts.ess_per_sec_min, " / ", result.nuts.ess_per_sec_mean)
        end
        if get(result, :advi, nothing) !== nothing
            println("  ADVI elapsed_sec: ", result.advi.elapsed_sec)
            haskey(result.advi, :gradient_used) && println("  ADVI gradient_used: ", result.advi.gradient_used)
            haskey(result.advi, :steady_state_solver_used) && println("  ADVI steady_state_solver_used: ", result.advi.steady_state_solver_used)
            haskey(result.advi, :zygote_fallback_used) && println("  ADVI zygote_fallback_used: ", result.advi.zygote_fallback_used)
            haskey(result.advi, :optimization_converged) && println("  ADVI optimization_converged: ", result.advi.optimization_converged)
            haskey(result.advi, :optimizer_iterations) && println("  ADVI optimizer_iterations: ", result.advi.optimizer_iterations)
            haskey(result.advi, :neg_elbo_min) && println("  ADVI neg_elbo_min: ", result.advi.neg_elbo_min)
            haskey(result.advi, :measures) && println("  ADVI WAIC: ", result.advi.measures.waic)
            haskey(result.advi, :posterior_samples) && println("  ADVI posterior_samples: ", result.advi.posterior_samples)
            haskey(result.advi, :posterior_mean) && println("  ADVI posterior_mean: ", result.advi.posterior_mean)
            haskey(result.advi, :posterior_std) && println("  ADVI posterior_std: ", result.advi.posterior_std)
            haskey(result.advi, :posterior_q) && println("  ADVI posterior_q (2.5%, 50%, 97.5%): ", result.advi.posterior_q)
        end
        if get(result, :mh, nothing) !== nothing && get(result, :nuts, nothing) !== nothing &&
           haskey(result.mh, :ess_per_sec_min) && haskey(result.nuts, :ess_per_sec_min)
            println("  MH/NUTS ESS-per-sec ratio (min): ", result.mh.ess_per_sec_min / result.nuts.ess_per_sec_min)
        end
    elseif get(result, :compare, nothing) === :trace_subset_gradient_three_way
        println("compare_trace_subset_gradient_benchmarks:")
        println("  time_forwarddiff_sec: ", result.time_forwarddiff_sec)
        println("  time_finitediff_sec: ", result.time_finitediff_sec)
        println("  time_zygote_sec: ", result.time_zygote_sec)
        println("  fastest_method: ", result.fastest_method)
        println("  forwarddiff_vs_finitediff max_abs / max_rel: ", result.forwarddiff_vs_finitediff_max_abs_diff, " / ", result.forwarddiff_vs_finitediff_max_rel_diff)
        println("  forwarddiff_vs_zygote max_abs / max_rel: ", result.forwarddiff_vs_zygote_max_abs_diff, " / ", result.forwarddiff_vs_zygote_max_rel_diff)
        println("  finitediff_vs_zygote max_abs / max_rel: ", result.finitediff_vs_zygote_max_abs_diff, " / ", result.finitediff_vs_zygote_max_rel_diff)
        haskey(result, :fd_ε) && result.fd_ε !== nothing && println("  fd_ε: ", result.fd_ε)
        haskey(result, :param_indices) && println("  param_indices: ", result.param_indices)
        haskey(result, :nparams) && println("  nparams: ", result.nparams)
        haskey(result, :θlen) && println("  θlen: ", result.θlen)
        haskey(result, :hmm_checkpoint_steps) && println("  hmm_checkpoint_steps: ", result.hmm_checkpoint_steps)
    elseif get(result, :compare, nothing) === :forwarddiff_vs_finitediff
        println("benchmark_trace_compare_forwarddiff_vs_finitediff:")
        println("  time_finitediff_sec: ", result.time_finitediff_sec)
        haskey(result, :gradient_check_reference) && println("  gradient_check_reference: ", result.gradient_check_reference)
        if result.time_forwarddiff_sec === nothing
            println("  time_forwarddiff_sec: (skipped) ", get(result, :forwarddiff_unavailable_reason, ""))
            println("  faster: finite-difference only")
        else
            println("  time_forwarddiff_sec: ", result.time_forwarddiff_sec)
            ρ = result.ratio_forwarddiff_over_finitediff
            println("  ratio (ForwardDiff time / finite-diff time): ", ρ)
            println("  faster: ", result.faster, "  (ratio < 1 ⇒ ForwardDiff quicker; ratio > 1 ⇒ finite-diff quicker)")
        end
        if haskey(result, :gradient_max_rel_diff) && result.gradient_max_rel_diff !== nothing
            ref = get(result, :gradient_check_reference, :forwarddiff_vs_fd)
            label = ref === :zygote_vs_fd ? "Zygote vs central FD" : "ForwardDiff vs central FD"
            println("  gradient_max_abs_diff ($label): ", result.gradient_max_abs_diff)
            println("  gradient_max_rel_diff: ", result.gradient_max_rel_diff)
        elseif get(result, :gradient_check_reference, nothing) === :unavailable
            println("  gradient agreement: skipped (coupled trace; pass check_zygote_vs_fd=true for Zygote vs FD)")
        elseif get(result, :gradient_check_reference, nothing) === :failed
            println("  gradient agreement: Zygote vs FD check failed (see @warn above)")
        end
        println("  param_indices: ", result.param_indices)
        println("  nparams: ", result.nparams)
        println("  n_ll_evals_per_grad: ", result.n_ll_evals_per_grad)
        println("  fd_ε: ", result.fd_ε)
        haskey(result, :θlen) && println("  θlen: ", result.θlen)
        haskey(result, :hmm_checkpoint_steps) && println("  hmm_checkpoint_steps: ", result.hmm_checkpoint_steps)
    elseif haskey(result, :elapsed_sec)
        println("elapsed_sec: ", result.elapsed_sec)
    elseif haskey(result, :time_no_checkpoint)
        println("benchmark_trace_zygote_gradient:")
        haskey(result, :success) && println("  success: ", result.success)
        haskey(result, :time_no_checkpoint) && println("  time_no_checkpoint: ", result.time_no_checkpoint)
        haskey(result, :time_with_checkpoint) && println("  time_with_checkpoint: ", result.time_with_checkpoint)
        haskey(result, :ratio) && println("  ratio (with / no checkpoint): ", result.ratio)
        haskey(result, :θlen) && println("  θlen: ", result.θlen)
        if haskey(result, :longest_timesteps)
            print("  longest_timesteps: ", result.longest_timesteps)
            if result.longest_timesteps === nothing
                print("  (not inferred — re-run `benchmark_trace_zygote_gradient` after updating StochasticGene, or trace layout has no `nframes`/trial length hint)")
            end
            println()
        end
        haskey(result, :checkpoint_chunks) && println("  checkpoint_chunks: ", result.checkpoint_chunks)
        haskey(result, :error) && println("  error: ", result.error)
    elseif haskey(result, :time_sec) && haskey(result, :param_indices)
        println("trace likelihood subset gradient (ForwardDiff or finite-diff):")
        println("  time_sec: ", result.time_sec)
        println("  param_indices: ", result.param_indices)
        println("  nparams: ", result.nparams)
        haskey(result, :θlen) && println("  θlen: ", result.θlen)
        haskey(result, :n_ll_evals_per_grad) && println("  n_ll_evals_per_grad: ", result.n_ll_evals_per_grad)
        haskey(result, :fd_ε) && println("  fd_ε: ", result.fd_ε)
        haskey(result, :hmm_checkpoint_steps) && println("  hmm_checkpoint_steps: ", result.hmm_checkpoint_steps)
    else
        println("(no standard benchmark fields; keys: ", join(string.(keys(result)), ", "), ")")
    end
    if haskey(result, :min_ess)
        println("min_ess: ", result.min_ess, "  mean_ess: ", result.mean_ess)
        println("max_rhat: ", result.max_rhat)
        println("ESS/s (min): ", result.ess_per_sec_min, "  ESS/s (mean): ", result.ess_per_sec_mean)
    end
    if haskey(result, :neg_elbo_min)
        println("ADVI Optim minimum (neg ELBO objective): ", result.neg_elbo_min)
    end
    return nothing
end
=#

#= Legacy NUTS/ADVI smoke tests moved to `test_features.jl`.

# ════════════════════════════════════════════════════════════════════════════════════
# SECTION: NUTS/ADVI Test Functions (for gradient-based Bayesian inference)
# ════════════════════════════════════════════════════════════════════════════════════

"""
    test_nuts_fit_smoke(; seed=42, n_samples=12, n_adapts=12, gradient=:finite)

Simple NUTS smoke test on a small synthetic RNA histogram problem. Validates that
`run_nuts_fit` returns proper structure (fits, stats, measures, nuts_info) matching
`run_mh` convention.

# Keywords
- `seed`: RNG seed
- `n_samples`: post-warmup samples
- `n_adapts`: adaptation steps
- `gradient`: gradient mode (`:finite`, `:ForwardDiff`, `:Zygote`)

# Returns
- `(fits, stats, measures, nuts_info)` from [`run_nuts_fit`](@ref)
"""
function test_nuts_fit_smoke(; seed::Int=42, n_samples::Int=12, n_adapts::Int=12, gradient::Symbol=:finite)
    gradient == :finite || @warn "test_nuts_fit_smoke delegates to test_run_nuts_fit_smoke, which uses gradient=:finite. Pass gradient=:finite for the existing smoke test path."
    return test_run_nuts_fit_smoke(; seed=seed, n_samples=n_samples, n_adapts=n_adapts)
end

"""
    test_run_advi_fit_smoke(; seed=42, maxiter=100, gradient=:finite)

Short ADVI smoke test for interactive validation. Reuses the existing synthetic
benchmark scenario from [`benchmark_inference_simrna_small`](@ref).
Returns `true` when the returned fit/statistics layout is consistent.
"""
function test_run_advi_fit_smoke(; seed::Int=42, maxiter::Int=100, gradient::Symbol=:finite)
    scen = benchmark_inference_simrna_small(seed=seed)
    opts = ADVIOptions(;
        maxiter=maxiter,
        n_mc=4,
        σ_floor=1e-4,
        init_s_raw=-4.0,
        verbose=false,
        gradient=gradient,
        time_limit=nothing,
    )
    fits, stats, measures, advi_info = run_advi_fit(
        scen.data, scen.model, opts;
        rng=MersenneTwister(seed),
        steady_state_solver=:augmented,
    )
    size(fits.param, 2) == 1 || return false
    length(fits.ll) == 1 || return false
    length(stats.meanparam) == size(fits.param, 1) || return false
    measures.waic isa Tuple || return false
    haskey(advi_info, :μ) || return false
    haskey(advi_info, :σ) || return false
    haskey(advi_info, :vb_elbo) || return false

    bench = benchmark_inference_run_advi(
        scen;
        rng=MersenneTwister(seed),
        advi_options=opts,
    )
    haskey(bench, :fits) || return false
    haskey(bench, :stats) || return false
    haskey(bench, :measures) || return false
    haskey(bench, :advi_info) || return false
    haskey(bench, :neg_elbo_min) || return false
    return true
end

"""
    test_run_advi_trace_smoke(; seed=42, maxiter=40, totaltime=30.0, ntrials=1)

Short trace-data ADVI smoke test using the existing GR2R2 benchmark scenario.
This exercises the AD-friendly trace likelihood path and the automatic switch
from `gradient=:Zygote` to the safer finite-difference route for traces.
"""
function test_run_advi_trace_smoke(; seed::Int=42, maxiter::Int=40, totaltime::Float64=30.0, ntrials::Int=1)
    scen = benchmark_inference_trace_gr2r2(; seed=seed, totaltime=totaltime, ntrials=ntrials)
    out = benchmark_inference_run_advi(
        scen;
        rng=MersenneTwister(seed),
        advi_options=ADVIOptions(
            ;
            maxiter=maxiter,
            n_mc=2,
            σ_floor=1e-4,
            init_s_raw=-4.0,
            verbose=false,
            gradient=:Zygote,
            time_limit=60.0,
        ),
    )
    isfinite(out.neg_elbo_min) || return false
    length(out.μ) == length(get_param(scen.model)) || return false
    length(out.σ) == length(get_param(scen.model)) || return false
    return true
end

"""
    test_nuts_vs_advi_consistency(; seed=42, n_samples=50)

Run both NUTS and ADVI on the same data/model and compare posterior estimates.
Both methods should produce similar posterior means (within MCMC variance).

# Returns
- `true` if consistency check passes
"""
function test_nuts_vs_advi_consistency(; seed::Int=42, n_samples::Int=50)
    scen = benchmark_inference_simrna_small(seed=seed)
    nuts = benchmark_inference_run_nuts_parallel(
        scen;
        nchains=1,
        n_samples=n_samples,
        n_adapts=max(20, n_samples ÷ 2),
        nuts_gradient=:finite,
        nuts_fd_ε=1e-4,
        verbose=false,
        progress=false,
        report_chains=false,
    )
    advi = benchmark_inference_run_advi(
        scen;
        rng=MersenneTwister(seed),
        advi_options=ADVIOptions(
            ;
            maxiter=200,
            n_mc=4,
            σ_floor=1e-4,
            init_s_raw=-4.0,
            verbose=false,
            gradient=:finite,
            time_limit=60.0,
        ),
    )
    mean_diff = norm(nuts.stats.meanparam .- inverse_transform_params(reshape(advi.μ, :, 1), scen.model)[:, 1])
    mean_diff < 2.0 || return false
    return true
end
=#
# you might also need to make `set_d`