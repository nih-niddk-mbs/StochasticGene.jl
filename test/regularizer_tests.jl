using Test, Random, LinearAlgebra, StochasticGene

@testset "Optional joint residence regularization" begin
    SG = StochasticGene
    covariance = (Matrix{Float64}(I, 4, 4) .- 0.25) .* 0.001
    target = [0.1, 0.2, 0.3, 0.4]
    term = JointResidenceRegularizer(target; covariance=covariance)
    @test_throws ArgumentError JointResidenceRegularizer([1, 1, 0, 0]; covariance=covariance)
    @test_throws ArgumentError JointResidenceRegularizer(target; covariance=-covariance)
    @test_throws ArgumentError JointResidenceRegularizer(target; covariance=zeros(4, 4))
    @test_throws ArgumentError JointResidenceRegularizer(target; covariance=covariance, units=(1, 1))
    @test_throws ArgumentError JointResidenceRegularizer(target; covariance=covariance, strength=-1)
    @test_throws ArgumentError SG.MHOptions(2, 0, 10., 1.; regularizers=("bad",))
    @test_throws ArgumentError SG.load_options(Dict(:inference_method=>:nuts, :regularizers=>(term,)))
    @test eval(Meta.parse(repr(term))).target == term.target

    rng = MersenneTwister(47)
    transitions = (([1, 2], [2, 1]), ([1, 2], [2, 1]))
    coupling = ((1, 2), [(2, 2, 1, 1)], [:free])
    rates = vcat([.1, .2, .3, .4, 1., 0., .2, 1., .3],
        [.15, .25, .4, .5, 1., 0., .2, 1., .3], [.5])
    traces = [rand(rng, 8, 2), rand(rng, 10, 2)]
    data = SG.TraceData("regularizer-test", "test", 1., (traces, Float64[], [0., 0.], [1., 1.], [0., 0.]))
    _, model = SG.make_trace_datamodel(data, rates, transitions, (2, 2), (1, 1), (0, 0), (1, 1),
        (SG.prob_Gaussian, SG.prob_Gaussian), 4, "", true, false, coupling)
    param = SG.get_param(model)
    raw = SG.loglikelihood(param, data, model)
    report = regularization_report(param, data, model; regularizers=(term,))
    @test report.loglikelihood == raw[1]
    @test report.pointwise == raw[2]
    @test report.logprior == SG.logprior(param, model)
    @test sum(only(report.fractions)) ≈ 1
    delta = only(report.fractions) .- target
    @test report.penalty ≈ dot(delta, pinv(covariance) * delta) / 2
    @test report.objective ≈ raw[1] + report.logprior - report.penalty
    doubled = regularization_report(param, data, model; regularizers=(term, term))
    @test doubled.penalty ≈ 2report.penalty
    SG._with_regularizers((term,), data, model) do
        proposed = param .+ 0.02
        current_prior = SG._mh_sampling_prior(param, data, model, SG.MHOptions(2, 0, 10., 1.; regularizers=(term,)))
        @test current_prior ≈ report.logprior - report.penalty
        proposed_report = regularization_report(proposed, data, model; regularizers=term)
        _, _, proposed_penalty = SG._mh_evaluate(proposed, data, model,
            SG.MHOptions(2, 0, 10., 1.; regularizers=(term,)))
        @test proposed_penalty ≈ proposed_report.penalty
        @test current_prior == SG._mh_sampling_prior(param, data, model,
            SG.MHOptions(2, 0, 10., 1.; regularizers=(term,)))
    end
    matching = JointResidenceRegularizer(only(report.fractions); covariance=covariance)
    @test regularization_report(param, data, model; regularizers=matching).penalty < 1e-20
    reverse_term = JointResidenceRegularizer(target; covariance=covariance, units=(2, 1))
    reversed = regularization_report(param, data, model; regularizers=reverse_term)
    @test only(reversed.fractions) ≈ only(report.fractions)[[1, 3, 2, 4]]
    @test !haskey(task_local_storage(), :stochasticgene_regularizers)
    tasks = [Threads.@spawn regularization_report(param, data, model; regularizers=(term,)) for _ in 1:2]
    @test fetch(tasks[1]).penalty == fetch(tasks[2]).penalty
    @test fetch(tasks[1]).fractions[1] !== fetch(tasks[2]).fractions[1]
    @test_throws ErrorException SG._with_regularizers((term,), data, model) do
        error("cleanup")
    end
    @test !haskey(task_local_storage(), :stochasticgene_regularizers)
    zero = JointResidenceRegularizer(target; covariance=covariance, strength=0)
    empty_options = SG.MHOptions(20, 5, 60., 1.)
    zero_options = SG.MHOptions(20, 5, 60., 1.; regularizers=(zero,))
    @test isempty(zero_options.regularizers)
    Random.seed!(92)
    unregularized, _ = SG.metropolis_hastings(data, model, empty_options)
    Random.seed!(92)
    disabled, _ = SG.metropolis_hastings(data, model, zero_options)
    @test unregularized.param == disabled.param
    @test unregularized.ll == disabled.ll
    @test unregularized.accept == disabled.accept
    options = SG.MHOptions(20, 5, 60., 1.; regularizers=(term,))
    Random.seed!(92)
    fit, _ = SG.metropolis_hastings(data, model, options)
    @test fit.total == 20
    @test fit.prior == SG.logprior(fit.param[:, end], model)
    @test fit.llml == SG.loglikelihood(fit.parml, data, model)[1]
    for i in axes(fit.param, 2)
        @test fit.ll[i] == SG.loglikelihood(fit.param[:, i], data, model)[1]
    end
    mktempdir() do folder
        spec = SG.fit_default_spec()
        spec[:regularizers] = (term,)
        file = SG.write_fitfile_full(joinpath(folder, "fit.jl"), spec)
        script = read(file, String)
        @test occursin("regularizers=", script)
        @test occursin("JointResidenceRegularizer", script)
        @test Meta.parseall(script) isa Expr
        SG.write_run_spec_jld2(joinpath(folder, "info_test.toml"), spec)
        loaded = SG.read_run_spec(joinpath(folder, "info_test.jld2"))
        @test only(loaded[:regularizers]).target == term.target
        @test only(loaded[:regularizers]).covariance == term.covariance
        SG._write_regularization_report(folder, data, model, fit,
            (medparam=SG.inverse_transform_params(fit.param[:, end], model),), options)
        @test length(filter(f -> startswith(f, "regularization_") && endswith(f, ".csv"), readdir(folder))) == 1
    end
end
