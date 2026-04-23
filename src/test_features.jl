# This file is part of StochasticGene.jl
#
# test_features.jl
#
# Feature-level correctness, compatibility, and regression checks.
# Loaded by `src/test.jl` so these `test_*` helpers remain available directly after
# `using StochasticGene`, while keeping the full-stack suite in `src/test.jl` leaner.

# ════════════════════════════════════════════════════════════════════════════════════
# TRACE / COUPLING FEATURE CORRECTNESS
# ════════════════════════════════════════════════════════════════════════════════════

"""
    test_trace_full(; coupling, G, R, S, insertstep, transitions, rtarget, noisepriors,
                    interval, totaltime, ntrials, method)

Simulate coupled joint traces then compare log-likelihood evaluated by the legacy RG-stack and the
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
(default: `test_compare_3unit`), prints per-histogram sums and a few bins so you can
see which histogram (ON/OFF/ONG/OFFG, which unit) is off, and returns
`(true, cme_vec, sim_vec)`.
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

# ════════════════════════════════════════════════════════════════════════════════════
# COMPATIBILITY / REGRESSION / AD FEATURE TESTS
# ════════════════════════════════════════════════════════════════════════════════════

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
        :maxtime        => 30.0,
        :label          => "tracejoint-HBEC-nstate_enhancer-gene11",
        :nchains        => 1,
        :samplesteps    => 100000,
        :datatype       => "tracejoint",
        :warmupsteps    => 0,
        :cell           => "HBEC",
        :resultfolder   => "3Prime-coupled-test",
        :propcv         => 0.05,
        :ratetype       => "median",
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
        :trace_specs    => default_trace_specs_for_coupled((1.6666666666666667, 1.0, -1.0), Bool[1, 1], 2),
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

# ════════════════════════════════════════════════════════════════════════════════════
# GRADIENT-BASED INFERENCE FEATURE SMOKE TESTS
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
