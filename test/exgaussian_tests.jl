using Test, StochasticGene, Random, Statistics

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
