# This file is part of StochasticGene.jl
#
# runtests.jl
#
# Role: **orchestrate** `Pkg.test()` and **assert** outcomes. Reusable scenarios and
# computations called from this file live in `src/test.jl` as exported `test_*` helpers,
# so each check can also be rerun interactively (e.g. `using StochasticGene; test_fit_simrna()`).
#
# Standalone files under `test/` are developer-focused and are not included here.
#
# Default `Pkg.test()` runs AD checks and other fast tests only. Long MCMC / fit
# integration tests are gated: set environment variable
#   STOCHASTICGENE_FULL_TESTS=1
# when validating the full stack (e.g. after end-to-end AD is confirmed).

using Random
using StochasticGene
using Test

const FULL_TESTS = get(ENV, "STOCHASTICGENE_FULL_TESTS", "0") == "1"

@testset "StochasticGene" begin

    @testset "maxtime_seconds" begin
        @test StochasticGene.maxtime_seconds(3600) == 3600.0
        @test StochasticGene.maxtime_seconds(3600.0) == 3600.0
        @test StochasticGene.maxtime_seconds("2h") == 7200.0
        @test StochasticGene.maxtime_seconds("90m") == 5400.0
        @test StochasticGene.maxtime_seconds("120") == 120.0
        @test StochasticGene.maxtime_seconds(" 1.5h ") == 5400.0
    end

    @testset "Biowulf script / swarm emission" begin
        mktempdir() do dir
            script = joinpath(dir, "gene_fit.jl")
            trace_row = (unit=1, interval=1.0, start=1.0, t_end=-1.0, zeromedian=true, background=0.0, active_fraction=1.0)
            StochasticGene.write_fitfile_genes(
                script, 2, (:rna, :rnadwelltime), (rna="data/rna/", dwelltime=("ON.csv", "OFF.csv")),
                "CELL", "MOCK", "res", "lb",
                Int[], tuple(), ([1, 2], [2, 1]), 2, 0, 0, 1, tuple(), nothing, ".", 60.0, 6.0,
                Float64[], 1, 10.0, Int[], -1.0, "", StochasticGene.prob_Gaussian, [], tuple(), "median", 0.01,
                1000, 0, 1.0, 1.0, false, false, false, "Tsit5()", "", true, 3, 1, 0.05,
                [trace_row], [];
                inference_method=:mh,
            )
            s = read(script, String)
            @test occursin("@time fit(;", s)
            @test occursin("gene=ARGS[1]", s)
            @test !occursin("fit(2,", s)
            @test !occursin("infolder", s) && !occursin("inlabel", s)
            @test occursin(":rna", s) && occursin(":rnadwelltime", s)
            swarmf = joinpath(dir, "batch.swarm")
            StochasticGene.emit_swarm_batch(swarmf, 2, 1, "gene_fit.jl", ["GENE1", "GENE2"])
            sw = read(swarmf, String)
            @test occursin("gene_fit.jl", sw) && occursin("GENE1", sw) && occursin("GENE2", sw)
            staged = joinpath(dir, "staged")
            StochasticGene.stage_run(joinpath(dir, "src"), staged; copy_specs=false)
            @test isdir(joinpath(staged, "specs"))
            @test isfile(joinpath(staged, "inputs", "manifest.toml"))
            mg = StochasticGene.model_grid(Gset=[2], Rset=[0], Sset=[0], insertset=[1])
            @test length(mg) == 1 && mg[1].G == 2
            outm = joinpath(dir, "minimal")
            mkpath(outm)
            StochasticGene.makeswarm(; filedir=outm, key="k1", juliafile="once.jl", resultfolder="rf", root=".", nchains=3)
            ms = read(joinpath(outm, "once.jl"), String)
            @test occursin("key=", ms) && occursin("\"k1\"", ms)
            sw1 = read(joinpath(outm, "fit.swarm"), String)
            @test occursin("once.jl", sw1)
            rb = joinpath(dir, "results", "batchres")
            mkpath(rb)
            touch(joinpath(rb, "info_mykey.jld2"))
            bout = joinpath(dir, "batchout")
            kfound = StochasticGene.makeswarm_folder("batchres"; root=dir, filedir=bout)
            @test kfound == ["mykey"] || sort(kfound) == ["mykey"]
            @test isfile(joinpath(bout, "fitscript_mykey.jl"))
            @test isfile(joinpath(bout, "fit.swarm"))
        end
    end

    @testset "RNA + ON/OFF + dwell" begin
        h1, h2 = StochasticGene.test_fit_rnaonoff()
        @test isapprox(h1, h2, rtol=0.05)

        h1, h2 = StochasticGene.test_fit_rnadwelltime()
        @test isapprox(h1, h2, rtol=0.3)

        @test StochasticGene.test_load_model_keyword_compatibility()
    end

    @testset "get_rates_ad consistency" begin
        @test StochasticGene.test_get_rates_ad_consistency()
    end

    @testset "run_nuts_fit smoke" begin
        @test StochasticGene.test_run_nuts_fit_smoke()
    end

    if FULL_TESTS
        @testset "RNA histograms" begin
            h1, h2 = StochasticGene.test_compare()
            @test isapprox(h1, h2, rtol=0.05)

            h1, h2 = StochasticGene.test_fit_simrna()
            @test isapprox(h1, h2, rtol=0.05)

            h1, h2 = StochasticGene.test_fit_rna()
            @test isapprox(h1, h2, rtol=0.05)
        end

        @testset "RNA + ON/OFF + dwell" begin
            h1, h2 = StochasticGene.test_fit_rnaonoff()
            @test isapprox(h1, h2, rtol=0.05)

            h1, h2 = StochasticGene.test_fit_rnadwelltime()
            @test isapprox(h1, h2, rtol=0.3)
        end

        @testset "Traces (single-unit)" begin
            lower, target, upper = StochasticGene.test_fit_trace()
            @test all(lower .<= target .<= upper)
            # `test_fit_trace_hierarchical` fits a very high-dimensional hierarchical
            # state (slow, machine-dependent sample counts under maxtime). Run it
            # interactively when validating hierarchical traces; see `src/test.jl`.
        end

        @testset "Coupled traces and dwell" begin
            lower, target, upper = StochasticGene.test_fit_tracejoint()
            @test all(lower .<= target .<= upper)

            # Basic coupled dwell-time vs simulator comparison (short 3-unit example)
            cme_vec, sim_vec = StochasticGene.test_compare_3unit()
            @test length(cme_vec) == length(sim_vec)
        end

        @testset "MH vs NUTS posterior (GM, simulated RNA histogram)" begin
            # Same likelihood and priors; compare posterior mean rates (Stats.meanparam).
            # Large MH chain + moderate NUTS; finite-diff gradients for robust NUTS.
            res = StochasticGene.test_fit_simrna_mh_nuts(;
                totalsteps=80_000,
                sim_seed=42,
                mh_samplesteps=35_000,
                mh_warmup=8_000,
                mh_maxtime=240.0,
                mh_seed=101,
                nuts_n_samples=1_200,
                nuts_n_adapts=600,
                nuts_gradient=:finite,
                nuts_fd_ε=1e-4,
                rng_nuts=MersenneTwister(202),
            )
            @test length(res.mean_mh) == length(res.mean_nuts)
            # Monte Carlo noise: allow generous margin on rate-space means
            @test res.max_abs_diff < 1.25
        end
    end

    @testset "trace_specs utilities" begin
        @test StochasticGene.test_trace_specs_utilities()
    end

    @testset "CombinedData (tracerna split / round-trip + canonical key order)" begin
        len = 5
        h = fill(0.1, len)
        trace = [rand(4, 2)]
        tup = (trace, 0.0, 0.0, 4, 1.0)
        d = StochasticGene.TraceRNAData("lab", "gene", 1.0, tup, len, h, 1.0, Int[])
        b = StochasticGene.CombinedData(d)
        @test StochasticGene.combined_modalities(b) == (:rna, :trace)
        td = StochasticGene.TraceData("lab", "gene", 1.0, tup, Int[])
        rna = StochasticGene.RNAData("lab", "gene", len, h, 1.0, Int[])
        b2 = StochasticGene.CombinedData((trace=td, rna=rna))
        @test StochasticGene.combined_modalities(b2) == (:rna, :trace)
        @test typeof(b) === typeof(b2)
        @test StochasticGene.normalize_datatype((:trace, :rna)) == (:rna, :trace)
        @test StochasticGene.normalize_datatype(["trace", "rna"]) == (:rna, :trace)
        @test_throws ArgumentError StochasticGene.normalize_datatype((:rna, :rna))
        d2 = StochasticGene.reconstruct_tracerna(b)
        @test d2.label == d.label && d2.gene == d.gene && d2.interval == d.interval
        @test d2.nRNA == d.nRNA && d2.histRNA == d.histRNA && d2.yield == d.yield
        @test d2.trace == d.trace && d2.units == d.units

        mktempdir() do root
            mkpath(joinpath(root, "data", "rna"))
            mkpath(joinpath(root, "data", "dwelltime"))
            touch(joinpath(root, "data", "dwelltime", "on.csv"))
            touch(joinpath(root, "data", "dwelltime", "off.csv"))
            raw = (
                rna = "rna",
                dwelltime = ["dwelltime/on.csv", "dwelltime/off.csv"],
            )
            resolved = StochasticGene._data_folder_path(raw, root)
            @test resolved.rna == joinpath(root, "data", "rna")
            @test resolved.dwelltime == [
                joinpath(root, "data", "dwelltime", "on.csv"),
                joinpath(root, "data", "dwelltime", "off.csv"),
            ]
            @test raw.dwelltime == ["dwelltime/on.csv", "dwelltime/off.csv"]
        end
    end

    @testset "CombinedData RNA+dwell reporter components" begin
        transitions = ([1, 2], [2, 1], [2, 3], [3, 2])
        G, R, S, insertstep = 3, 2, 2, 1
        onstates = [Int[], Int[], [2, 3], [2, 3]]
        bins = [collect(0.1:0.1:2.0) for _ in 1:4]
        dwell = [ones(length(bins[1])) for _ in 1:4]
        dwell_data = StochasticGene.DwellTimeData("test", "test", bins, dwell, ["ON", "OFF", "ONG", "OFFG"])
        rna_data = StochasticGene.RNAData("test", "test", 20, ones(20), 1.0, Int[])
        combined = StochasticGene.CombinedData((rna=rna_data, dwelltime=dwell_data))

        reporter, components = StochasticGene.make_reporter_components(
            combined, transitions, G, R, S, insertstep, "", onstates, 0.01,
            StochasticGene.prob_Gaussian, Float64[], tuple(), 1,
        )
        @test StochasticGene.combined_modalities(combined) == (:dwelltime, :rna)
        @test reporter !== nothing
        @test components isa StochasticGene.MTComponents

        r = [0.038, 0.3, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.00231]
        rmean = StochasticGene.prior_ratemean(
            transitions, R, S, insertstep, r[end], Float64[],
            StochasticGene.mean_elongationtime(r, transitions, R),
        )
        model = StochasticGene.load_model(
            combined, r, rmean, collect(1:length(r)-1), tuple(), transitions, G, R, S,
            insertstep, "", 2, 10.0, onstates, r[end], 0.05, StochasticGene.prob_Gaussian,
            Float64[], StochasticGene.Tsit5(), tuple(), tuple(), nothing,
        )
        θ = StochasticGene.get_param(model)
        ll, pred = StochasticGene.loglikelihood(θ, combined, model, StochasticGene.MHOptions(1, 0, 1.0, 1.0))
        @test isfinite(ll)
        @test length(pred) == length(StochasticGene.datapdf(combined))
        @test StochasticGene.filename(combined, model) == "_test_test_3221_2.txt"
    end

    @testset "CombinedData independent likelihood assembly" begin
        coupling = tuple()
        transitions = ([1, 2], [2, 1])
        G = 2
        R = 2
        S = 1
        insertstep = 1
        rtarget = Float64[
            0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 1.0, 0.0, 0.05, 1.0, 0.05,
        ]
        interval = 1.0
        # Fixed trace keeps this assembly test deterministic; short stochastic
        # simulations can occasionally produce empty traces.
        trace = [Float64[0.2, 1.1, 0.9, 0.3]]
        nRNA = 8
        h = fill(1 / nRNA, nRNA)
        data = StochasticGene.TraceRNAData(
            "tracerna", "test", interval, (trace, 0.0, 0.0, 4, 1.0),
            nRNA, h, 1.0, Int[],
        )
        priormean = StochasticGene.set_priormean(
            [], transitions, R, S, insertstep, 1.0, [0.0, 0.1, 1.0, 0.1],
            StochasticGene.mean_elongationtime(rtarget, transitions, R), tuple(), coupling, nothing,
        )
        model = StochasticGene.load_model(
            data, rtarget, priormean, Int[1], tuple(), transitions, G, R, S, insertstep,
            "", 1, 10.0, Int[], rtarget[StochasticGene.num_rates(transitions, R, S, insertstep)],
            0.2, StochasticGene.prob_Gaussian, [0.0, 0.1, 1.0, 0.1],
            StochasticGene.Tsit5(), tuple(), coupling, nothing,
        )
        θ = StochasticGene.get_param(model)
        combined = StochasticGene.CombinedData(data)
        opts = StochasticGene.MHOptions(1, 0, 1.0, 1.0)

        ll, pred = StochasticGene.loglikelihood(θ, combined, model, opts; steady_state_solver=:augmented)
        ll_rna, pred_rna = StochasticGene.loglikelihood(θ, StochasticGene.CombinedData((rna=combined.legs.rna,)), model, opts; steady_state_solver=:augmented)
        ll_trace, pred_trace = StochasticGene.loglikelihood(θ, StochasticGene.CombinedData((trace=combined.legs.trace,)), model, opts; steady_state_solver=:augmented)
        @test ll ≈ ll_rna + ll_trace
        @test pred ≈ vcat(pred_rna, pred_trace)

        ll_legacy, pred_legacy = StochasticGene.loglikelihood(θ, data, model; steady_state_solver=:augmented, hmm_stack=opts.likelihood_executor)
        @test ll ≈ ll_legacy
        @test pred ≈ pred_legacy

        ad_opts = StochasticGene.NUTSOptions(n_samples=1, n_adapts=1)
        ll_ad, pred_ad = StochasticGene.loglikelihood_ad(θ, combined, model, ad_opts; steady_state_solver=:augmented)
        ll_rna_ad, pred_rna_ad = StochasticGene.loglikelihood_ad(θ, StochasticGene.CombinedData((rna=combined.legs.rna,)), model, ad_opts; steady_state_solver=:augmented)
        ll_trace_ad, pred_trace_ad = StochasticGene.loglikelihood_ad(θ, StochasticGene.CombinedData((trace=combined.legs.trace,)), model, ad_opts; steady_state_solver=:augmented)
        @test ll_ad ≈ ll_rna_ad + ll_trace_ad
        @test pred_ad ≈ vcat(pred_rna_ad, pred_trace_ad)
    end

    @testset "inference benchmark (smoke)" begin
        scen = StochasticGene.benchmark_inference_simrna_small(seed=7)
        @test hasproperty(scen, :data) && hasproperty(scen, :model)
        @test length(StochasticGene.get_param(scen.model)) > 0
        tr = StochasticGene.benchmark_inference_trace_gr2r2(seed=11, totaltime=120.0, ntrials=2)
        @test hasproperty(tr, :data) && hasproperty(tr, :model)
        @test length(StochasticGene.get_param(tr.model)) > 0
        cj = StochasticGene.benchmark_inference_trace_coupled_3x3(seed=11, totaltime=120.0, ntrials=2)
        @test hasproperty(cj, :data) && hasproperty(cj, :model)
        @test length(StochasticGene.get_param(cj.model)) > 0
        cjh = StochasticGene.benchmark_inference_trace_coupled_3x3_g3r0(seed=11, totaltime=120.0, ntrials=2)
        @test hasproperty(cjh, :data) && hasproperty(cjh, :model)
        @test length(StochasticGene.get_param(cjh.model)) > 0
    end

    @testset "staging_key_segment (Result_Fields only)" begin
        r = Result_Fields("rates", "scRNA", "", "CANX", "3331", "2")
        @test staging_key_segment(r, :gene) == "CANX"
    end

    @testset "stage_label_to_key key_mode=:fields (Summary vs Result)" begin
        mktempdir() do dir
            src = joinpath(dir, "in")
            mkpath(src)
            write(joinpath(src, "measures_scRNA__3331.csv"), "x\n")
            write(joinpath(src, "rates_scRNA_CANX_3331_2.txt"), "y\n")
            dst = joinpath(dir, "out")
            r = stage_label_to_key(
                src, dst;
                root=dir,
                root_dst=dir,
                key_mode=:fields,
                key_fields=(:label, :model),
                dry_run=true,
            )
            @test r.n == 2
            km = Dict(basename(row.dst) => row.key for row in r.rows)
            @test km["measures_scRNA-3331.csv"] == "scRNA-3331"
            @test km["rates_scRNA-3331.txt"] == "scRNA-3331"
            # (:label, :gene): aggregate CSV uses third stem token only for tuple compatibility with per-run files — not a gene.
            r2 = stage_label_to_key(
                src, dst;
                root=dir,
                root_dst=dir,
                key_mode=:fields,
                key_fields=(:label, :gene),
                dry_run=true,
            )
            @test r2.n == 2
            km2 = Dict(basename(row.dst) => row.key for row in r2.rows)
            @test km2["measures_scRNA.csv"] == "scRNA"
            @test km2["rates_scRNA-CANX.txt"] == "scRNA-CANX"
            # Default keyword `key_fields` for `:fields` (discriminating per-gene + batch CSV friendly).
            rdef = stage_label_to_key(
                src, dst;
                root=dir,
                root_dst=dir,
                key_mode=:fields,
                dry_run=true,
            )
            @test rdef.n == 2
            kmdef = Dict(basename(row.dst) => row.key for row in rdef.rows)
            @test kmdef["measures_scRNA-3331.csv"] == "scRNA-3331"
            @test kmdef["rates_scRNA-CANX-3331-2.txt"] == "scRNA-CANX-3331-2"
        end
    end

    @testset "stage_combine_rates (number_of_parameters, new_params, write_out, force)" begin
        mktempdir() do dir
            p1 = joinpath(dir, "a.txt")
            p2 = joinpath(dir, "b.txt")
            d1 = [1.0 2.0 3.0 4.0; 10.0 20.0 30.0 40.0]
            h1 = ["u1_a1", "u1_a2", "u1_a3", "u1_a4"]
            d2 = [5.0 6.0; 50.0 60.0]
            h2 = ["u2_b1", "u2_b2"]
            write_rates_table(p1, d1, h1)
            write_rates_table(p2, d2, h2)
            outpath = joinpath(dir, "merged_custom.txt")
            r = stage_combine_rates(
                [p1, p2],
                outpath,
                [2, 1],
                [9.0, 0.5, 0.6],
                ["Hidden_1", "Coupling_1", "Coupling_2"],
            )
            @test r.nsets == 2
            @test r.dst == abspath(outpath)
            dout, hout = read_rates_table(r.dst)
            @test size(dout) == (2, 12)
            @test dout[1, :] == [1.0, 2.0, 5.0, 9.0, 0.5, 0.6, 3.0, 4.0, 6.0, 9.0, 0.5, 0.6]
            @test hout == [
                "u1_a1", "u1_a2", "u2_b1", "Hidden_1_1", "Coupling_1_1", "Coupling_2_1",
                "u1_a3", "u1_a4", "u2_b2", "Hidden_1_2", "Coupling_1_2", "Coupling_2_2",
            ]

            out2 = joinpath(dir, "no_extra.txt")
            stage_combine_rates([p1, p2], out2, [2, 1], Float64[])
            d2only, h2only = read_rates_table(out2)
            @test size(d2only) == (2, 6)
            @test h2only == ["u1_a1", "u1_a2", "u2_b1", "u1_a3", "u1_a4", "u2_b2"]

            out3 = joinpath(dir, "dry.txt")
            r3 = stage_combine_rates([p1, p2], out3, [2, 1], Float64[], String[], false, false)
            @test !isfile(out3)
            @test r3.dst == abspath(out3)

            outhier = joinpath(dir, "merged_hier.txt")
            stage_combine_rates(
                [p1, p2],
                outhier,
                [2, 1],
                [9.0, 0.5, 0.6],
                ["Hidden_1", "Coupling_1", "Coupling_2"],
                true,
            )
            dh, hh = read_rates_table(outhier)
            @test size(dh) == (2, 9)
            @test dh[1, :] == [1.0, 2.0, 5.0, 3.0, 4.0, 6.0, 9.0, 0.5, 0.6]
            @test hh == [
                "u1_a1", "u1_a2", "u2_b1", "u1_a3", "u1_a4", "u2_b2",
                "Hidden_1", "Coupling_1", "Coupling_2",
            ]
        end
    end

    @testset "stage_combine_rates_from_csv convenience driver" begin
        mktempdir() do dir
            p1 = joinpath(dir, "u1_rates.txt")
            p2 = joinpath(dir, "u2_rates.txt")
            write_rates_table(p1, [1.0 2.0; 10.0 20.0], ["u1_a1", "u1_a2"])
            write_rates_table(p2, [5.0; 50.0;;], ["u2_b1"])

            csv_path = joinpath(dir, "couplings.csv")
            df = DataFrame(
                Model_name=["k_free", "k_act", "k_inh"],
                c12=["free", ">0", "<0"],
            )
            CSV.write(csv_path, df)

            specs = stage_combine_rates_specs_from_csv(
                csv_path,
                ["c12"],
                ["Coupling_1"];
                coupling_mode_values=Dict(:free => 0.0, :activate => 0.2, :inhibit => -0.2),
                base_new_params=[9.0],
                base_new_labels=["Hidden_1"],
            )
            @test length(specs) == 3
            @test specs[1].new_params == [9.0, 0.0]
            @test specs[2].new_params == [9.0, 0.2]
            @test specs[3].new_params == [9.0, -0.2]

            out = joinpath(dir, "out")
            res = stage_combine_rates_from_csv(
                csv_path,
                [p1, p2],
                [2, 1],
                out;
                coupling_mode_values=Dict(:free => 0.0, :activate => 0.2, :inhibit => -0.2),
                base_new_params=[9.0],
                base_new_labels=["Hidden_1"],
            )
            @test res.n == 3
            d_act, h_act = read_rates_table(joinpath(out, "rates_k_act.txt"))
            @test d_act[1, :] == [1.0, 2.0, 5.0, 9.0, 0.2]
            @test h_act == ["u1_a1", "u1_a2", "u2_b1", "Hidden_1", "Coupling_1"]

            out_flag = joinpath(dir, "out_flag")
            res_flag = stage_combine_rates_from_csv(
                csv_path,
                [p1, p2],
                [2, 1],
                out_flag;
                key_flag="pilotA",
                coupling_mode_values=Dict(:free => 0.0, :activate => 0.2, :inhibit => -0.2),
                base_new_params=[9.0],
                base_new_labels=["Hidden_1"],
            )
            @test res_flag.n == 3
            @test all(startswith(row.key, "pilotA-") for row in res_flag.rows)
            @test isfile(joinpath(out_flag, "rates_pilotA-k_act.txt"))
        end
    end

    @testset "stage-native make_fitscript helpers" begin
        mktempdir() do dir
            csv_path = joinpath(dir, "keys.csv")
            CSV.write(csv_path, DataFrame(Model_name=["sA", "sB"]))

            scripts = make_fitscripts_from_csv(
                csv_path;
                filedir=dir,
                juliafile="fitscript",
                resultfolder="rf",
                root=".",
            )
            @test length(scripts) == 2
            @test all(isfile, scripts)

            one = make_fitscript(
                "single-key";
                filedir=dir,
                juliafile="fitscript",
                resultfolder="rf_single",
            )
            @test isfile(one)
            @test basename(one) == "fitscript_single-key.jl"

            cmd = build_julia_script_command(
                "fitscript_sA.jl";
                nthreads=4,
                nprocs=3,
                project=".",
                extra_flags=["--check-bounds=no"],
            )
            @test occursin("--project=.", cmd)
            @test occursin("-t 4", cmd)
            @test occursin("-p 3", cmd)
            @test occursin("fitscript_sA.jl", cmd)

            cfile = make_commandfile_from_csv(
                csv_path;
                filedir=dir,
                commandfile="fit.commands",
                juliafile="fitscript",
                nthreads=2,
                nprocs=5,
            )
            @test isfile(cfile)
            ctext = read(cfile, String)
            @test occursin("-t 2", ctext)
            @test occursin("-p 5", ctext)
            @test occursin("fitscript_sA.jl", ctext)
            @test occursin("fitscript_sB.jl", ctext)

            cfile_direct = make_commandfile(
                ["fitscript_sA.jl", "fitscript_sB.jl"];
                filedir=dir,
                commandfile="direct.commands",
                nthreads=1,
                nprocs=2,
            )
            @test isfile(cfile_direct)
            @test occursin("fitscript_sA.jl", read(cfile_direct, String))

            swarm = make_swarmfile_from_csv(
                csv_path;
                filedir=dir,
                swarmfile="stagefit",
                juliafile="fitscript",
            )
            @test isfile(swarm)
            stext = read(swarm, String)
            @test occursin("fitscript_sA.jl", stext)
            @test occursin("fitscript_sB.jl", stext)

            both = make_fitscripts_and_swarm_from_csv(
                csv_path;
                filedir=joinpath(dir, "both"),
                swarmfile="combo",
                juliafile="emit",
                resultfolder="rf2",
            )
            @test isfile(both.swarm)
            @test length(both.scripts) == 2

            both_cmd = make_fitscripts_and_commandfile_from_csv(
                csv_path;
                filedir=joinpath(dir, "both_cmd"),
                commandfile="combo.commands",
                juliafile="emit2",
                resultfolder="rf3",
            )
            @test isfile(both_cmd.commandfile)
            @test length(both_cmd.scripts) == 2
        end
    end

end
