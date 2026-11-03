using Test, StochasticGene, Random, Statistics, CSV, DataFrames

@testset "Centered exGaussian emissions" begin
    SG=StochasticGene
    p=[-.05,.1,.7,.2,.15]
    for n in 0:3
        d=prob_ExGaussian(p,n)
        @test mean(d) ≈ p[1]+n*p[3]
        @test var(d) ≈ p[2]^2+n*p[4]^2+p[5]^2
        flat=copy(p);flat[5]=0.
        for x in (-10.,-2.,0.,1.,20.)
            normal=SG.prob_Gaussian(flat[1:4],n)
            @test SG.logpdf(prob_ExGaussian(flat,n),x) ≈ SG.logpdf(normal,x)
            @test SG.cdf(prob_ExGaussian(flat,n),x) ≈ SG.cdf(normal,x)
        end
    end
    # Independent scipy.stats.exponnorm reference: K=1.5, loc=-.3, scale=.2.
    d=SG.ExGaussianNoise(-.3,.2,.3)
    xs=[-10.,-2.,0.,1.,5.,20.]
    logs=[-1179.7355948632762,-38.06710538711433,.20013675535745512,
        -2.9071383095017183,-16.240471640118507,-66.2404716401185]
    cdfs=[0.,6.73357158313211e-19,.5667218578452918,.9836104462141494,.9999999734555313,1.]
    @test SG.logpdf.(Ref(d),xs) ≈ logs atol=1e-10
    @test SG.cdf.(Ref(d),xs) ≈ cdfs atol=1e-13
    @test SG.cdf(d,-Inf)==0
    @test SG.cdf(d,Inf)==1
    @test SG.logpdf(d,Inf)==-Inf
    tiny=SG.ExGaussianNoise(0.,1.,1e-12)
    for x in (-20.,0.,20.)
        @test SG.logpdf(tiny,x) ≈ SG.logpdf(SG.Normal(),x) atol=1e-8
    end
    draws=rand(MersenneTwister(83),d,100_000)
    @test abs(mean(draws)-mean(d))<.006
    @test abs(var(draws)-var(d))<.006
    @test_throws DomainError prob_ExGaussian([p[1:4];-.1],0)
    @test_throws ArgumentError prob_ExGaussian(p[1:4],0)
    @test length(prob_ExGaussian(p,[0,1,2],2))==2
    reporter=SG.HMMReporter(5,Int[],prob_ExGaussian,0,Int[],collect(10:14))
    f=Function[];inverse=Function[];sigma=Function[]
    SG.rate_transforms!(f,inverse,sigma,9,reporter,true)
    @test f[10]===identity
    @test f[12]===identity
    @test f[14]===log
    @test inverse[14](f[14](p[5])) ≈ p[5]
    a=[.9 .1;.2 .8];p0=[2/3,1/3];ds=prob_ExGaussian(p,[0,1])
    traces=[[-.2,.1,.8,1.,.1]]
    old=SG._ll_hmm(a,p0,ds,traces;hmm_stack=:fast)
    new=SG._ll_hmm_batched(a,p0,ds,traces;cache=SG.BatchedHMMCache(),batchsize=4)
    @test old[1] ≈ new[1]
    @test old[2] ≈ new[2]
end

@testset "Emission counts across coupled analysis and simulation" begin
    SG = StochasticGene
    transitions = (([1,2], [2,1]), ([1,2], [2,1]), ([1,2], [2,1]))
    G, R, S, steps = (2,2,2), (1,1,0), (0,0,0), (1,1,0)
    coupling = ((1,2,3), [(3,2,1,1)], [:free])
    kinetics = [.1,.2,.3,.4,1.]
    gaussian = [0.,.2,1.,.3]
    exgaussian = [gaussian; .1]
    hidden = [.1,.2,.3,1.]
    rates = vcat(kinetics, exgaussian, kinetics, exgaussian, hidden, [.2])
    args = (rates, transitions, G, R, S, steps, SG.prob_ExGaussian, coupling, [0.,1.,2.])
    inferred = SG.correlation_functions(args...)
    explicit = SG.correlation_functions(args...; noise_per_unit=[5,5,0])
    reference = SG._correlation_functions_reference(args...)
    @test all(isapprox(a,b; atol=1e-10) for (a,b) in zip(inferred, explicit))
    @test all(isapprox(a,b; atol=1e-10) for (a,b) in zip(inferred, reference))
    @test SG._emission_noise_count(SG.prob_GaussianMixture) == 5
    @test SG._emission_noise_count(SG.prob_GaussianMixture_6) == 7

    functions = (SG.prob_Gaussian, SG.prob_ExGaussian, SG.prob_Gaussian)
    mixed_rates = vcat(kinetics, gaussian, kinetics, exgaussian, hidden, [.2])
    ctx = SG.build_correlation_context(mixed_rates, transitions, G, R, S, steps,
        functions, coupling, [0.,1.,2.])
    @test length.(ctx.noiseparams) == [4,5,0]
    vector_ctx = SG.build_correlation_context(mixed_rates, transitions, G, R, S, steps,
        collect(functions), coupling, [0.,1.,2.])
    @test vector_ctx.p0 ≈ ctx.p0
    spec = Dict{Symbol,Any}(:transitions=>transitions, :G=>G, :R=>R, :S=>S,
        :insertstep=>steps, :probfn=>functions, :coupling=>coupling, :key=>"mixed-h",
        :noisepriors=>(gaussian, exgaussian, Float64[]), :trace_specs=>[(interval=1.,)])
    @test SG._correlation_noise_counts(spec) == [4,5,0]
    @test SG._correlation_noise_counts(merge(spec, Dict(:noisepriors=>(gaussian, exgaussian)))) == [4,5,0]
    @test SG._correlation_noise_counts(merge(spec, Dict(:noisepriors=>exgaussian))) == [5,5,0]
    mktempdir() do dir
        SG.write_run_spec_jld2(joinpath(dir,"info_mixed-h.toml"), spec)
        file = joinpath(dir,"rates_mixed-h.txt")
        CSV.write(file, DataFrame(repeat(permutedims(mixed_rates),4,1), :auto))
        @test SG.write_correlation_functions_key(dir; lags=[0.,1.,2.], validate=true) == 1
        raw = SG.correlation_functions_file(file)
        @test all(isfinite, raw[9])
        @test SG.write_correlation_functions_key(dir; lags=[0.,1.,2.],
            trace_center=true, window_lengths=4) == 1
        outputs = filter(f -> occursin("trace-centered",f), readdir(dir; join=true))
        @test length(outputs) == 2
        @test all(all(isfinite, CSV.read(f,DataFrame).tau) for f in outputs)
    end
    Random.seed!(87)
    traces = SG.simulate_trace_vector(mixed_rates, transitions, G, R, S, steps,
        coupling, 1., 10., 2; probfn=functions, col=[2,3], warmuptime=0.)
    @test length(traces) == 2
    @test all(t -> size(t,2) == 4 && all(isfinite,t), traces)
    trial = SG.simulate_trials(mixed_rates, transitions, G, R, S, steps, coupling,
        2, 10., [0.,1.,2.]; probfn=functions, warmuptime=0.)
    @test all(isfinite, trial.ccON_theory)
    @test all(isfinite, trial.ccON_mean)
end

@testset "Mixed exGaussian compact cache" begin
    SG = StochasticGene
    rng = MersenneTwister(129)
    a = rand(rng,4,4)
    a ./= sum(a; dims=2)
    p0 = fill(.25,4)
    reporters = [SG.HMMReporter(4,[0,0,1,1],SG.prob_Gaussian,0,Int[],Int[]),
        SG.HMMReporter(5,[0,1,0,1],SG.prob_ExGaussian,0,Int[],Int[])]
    noise = [[[0.,.2,1.,.3], [0.,.2,1.,.3,.1]] for _ in 1:3]
    traces = [rand(rng,3+i,2) for i in 1:3]
    cache = SG.BatchedHMMCache()
    function check()
        reference = SG._ll_hmm(noise,a,p0,reporters,traces; hmm_stack=:fast)
        result = SG._ll_hmm_batched(noise,a,p0,reporters,traces; cache=cache)
        @test result[1] ≈ reference[1]
        @test result[2] ≈ reference[2]
    end
    for batchsize in (1,2,8)
        cache.batchsize = batchsize
        check()
    end
    @test length(cache.emissions) == 3
    stored = copy(cache.emissions)
    check()
    @test all(cache.emissions[i] === stored[i] for i in eachindex(stored))
    noise[2][2][5] *= 2
    check()
    @test cache.emissions[1] === stored[1]
    @test cache.emissions[2] !== stored[2]
    noise[2][2][5] = 0.
    check()
end
