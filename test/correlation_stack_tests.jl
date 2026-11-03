using Test, StochasticGene, CSV, DataFrames, LinearAlgebra

@testset "Unified correlation engine" begin
    SG = StochasticGene
    p = SG.get_3unit_model_params()
    lags = collect(0.0:5.0:30.0)
    args = (p.r, p.transitions, p.G, p.R, p.S, p.insertstep,
        SG.prob_Gaussian, p.coupling, lags)
    options = (observed_units=p.observed_units, noise_per_unit=p.noiseparams)
    ctx = SG.build_correlation_context(args...; options...)
    raw = SG.correlation_functions(args...; options...)
    reference = SG._correlation_functions_reference(args...; options...)
    @test all(isapprox(a, b; atol=1e-10, rtol=1e-8) for (a,b) in zip(raw, reference))
    pairs = [((:ON, 1), (:ON, 2)), ((:reporters, 1), (:reporters, 2)),
        ((:gene_state, 1, 2), (:gene_state, 2, 2))]
    serial = SG.correlate_observables(ctx, pairs)
    cached_count = length(ctx.propagated_observables)
    threaded = SG.correlate_observables(ctx, pairs; threaded=true, max_concurrent_jobs=2)
    @test serial == threaded
    @test length(ctx.propagated_observables) == cached_count
    custom_x = (kind=:custom, label="same", values=zeros(ctx.components.N))
    custom_y = (kind=:custom, label="same", values=ones(ctx.components.N))
    @test all(iszero, SG.correlate_observables(ctx, custom_x, custom_y).cc)
    @test SG.correlate_observables(ctx, custom_y, custom_y).cc ≈ ones(length(ctx.tau))
    # Returned arrays must not expose mutable cache entries.
    serial[1].cc[1] = -100.0
    @test SG.correlate_observables(ctx, pairs[1]...).cc == threaded[1].cc
    centered = SG.correlate_observables(ctx, pairs[1]...;
        trace_center=true, window_lengths=7)
    @test centered.cc == threaded[1].cc
    @test centered.cc_centered ≈ SG.finite_trace_center_correlation(centered.cc, 7)
    @test_throws ArgumentError SG.correlate_observables(ctx, pairs[1]...; trace_center=true)
    @test_throws ArgumentError SG._correlation_lag_steps([0.0, 1.0, 1.0])
    @test_throws ArgumentError SG._run_correlation_jobs(identity, [1]; max_concurrent_jobs=0)
    @test SG._run_correlation_jobs(x -> x^2, 1:12; max_concurrent_jobs=2) == (1:12).^2
    @test SG._correlation_default_pairs(nothing, nothing, [2,3]) ==
        [((:ON, 2), (:ON, 3)), ((:reporters, 2), (:reporters, 3))]
    # Check a direct free-coupling architecture independently of the hidden case.
    direct_args = (vcat(p.r[1:18], [-0.3, 2.0]), p.transitions[1:2], p.G[1:2],
        p.R[1:2], p.S[1:2], p.insertstep[1:2], SG.prob_Gaussian,
        ((1,2), [(1,2,2,1), (2,1,1,2)], [:free, :free]), lags)
    direct = SG.correlation_functions(direct_args...)
    direct_reference = SG._correlation_functions_reference(direct_args...)
    @test all(isapprox(a, b; atol=1e-10, rtol=1e-8) for (a,b) in zip(direct, direct_reference))

    mktempdir() do dir
        df = SG._correlation_general_dataframe(threaded)
        input = joinpath(dir, "crosscorrelation-general_test.csv")
        CSV.write(input, df)
        output = SG.write_correlation_functions_centered(input; window_lengths=7,
            window_interval=5.0)
        saved = CSV.read(output, DataFrame)
        @test names(saved) == names(df)
        @test saved.cc == df.cc
        @test saved.cc_centered[1:13] ≈ centered.cc_centered
        serial_paths = SG.write_correlation_functions_centered_folder(dir;
            window_lengths=7, window_interval=5.0, threaded=false)
        parallel_paths = SG.write_correlation_functions_centered_folder(dir;
            window_lengths=7, window_interval=5.0, threaded=true, max_concurrent_jobs=2)
        @test serial_paths == parallel_paths == [output]
    end

    @testset "File and folder frontends" begin
        mktempdir() do dir
            spec = Dict{Symbol,Any}(:transitions=>p.transitions, :G=>p.G, :R=>p.R,
                :S=>p.S, :insertstep=>p.insertstep, :probfn=>SG.prob_Gaussian,
                :coupling=>p.coupling, :trace_specs=>[(interval=5/3,)], :key=>"saved-key")
            for key in ("first", "second")
                SG.write_run_spec_jld2(joinpath(dir, "info_$(key).toml"), spec)
                CSV.write(joinpath(dir, "rates_$(key).txt"),
                    DataFrame(repeat(permutedims(p.r), 4, 1), :auto))
            end
            requested = collect(0:5/3:5)
            @test SG.write_correlation_functions_key(dir; lags=requested, trace_center=true,
                window_lengths=7, threaded=false, validate=true) == 2
            files = sort(filter(f -> endswith(f, ".csv"), readdir(dir; join=true)))
            @test length(files) == 4
            @test all(f -> occursin("trace-centered", f), files)
            serial_outputs = [CSV.read(f, DataFrame) for f in files]
            @test SG.write_correlation_functions_key(dir; lags=requested, trace_center=true,
                window_lengths=7, threaded=true, max_concurrent_jobs=2) == 2
            @test all(isequal(CSV.read(f, DataFrame), df) for (f, df) in zip(files, serial_outputs))
            # File adapter infers the same frame interval and uses the same engine.
            SG.write_correlation_functions_file(joinpath(dir, "rates_first.txt");
                lags=requested, trace_center=true, window_lengths=7)
            @test all(isequal(CSV.read(f, DataFrame), df) for (f, df) in zip(files, serial_outputs))
            @test SG.write_correlation_functions(dir; lags=requested, trace_center=true,
                window_lengths=7, max_concurrent_jobs=1) === nothing
        end
    end

    @testset "Empirical frame interval" begin
        mktempdir() do dir
            a = Float64[0, 0, 1, 2, 1, 0, 0, 1, 0, 0, 2, 1]
            b = circshift(a, 1)
            f1, f2, out = joinpath.(dir, ["unit1.csv", "unit2.csv", "empirical.csv"])
            CSV.write(f1, DataFrame(Reporters_1=a))
            CSV.write(f2, DataFrame(Reporters_1=b))
            dt = 5/3
            SG.write_correlation_functions_empirical(dir, f1, f2, out;
                lags=[0.0, dt, 2dt], traceinfo=(dt, 1.0, -1))
            saved = CSV.read(out, DataFrame)
            expected = SG.compute_correlation_functions([hcat(a,b)], [-2dt, -dt, 0.0, dt, 2dt];
                frame_interval=dt)
            @test saved.cc_Reporters ≈ expected.cc
        end
    end
end
