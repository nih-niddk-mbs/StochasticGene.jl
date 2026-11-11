# Included by runtests.jl after using StochasticGene.
using Test, SparseArrays, LinearAlgebra, ForwardDiff
using Zygote
const SGFPT = StochasticGene
struct FPTRegressionModel <: SGFPT.AbstractGeneTransitionModel end
SGFPT.predictedfn(param, data::SGFPT.DwellTimeData, ::FPTRegressionModel; kwargs...) = param
SGFPT.predictedfn(param, data::SGFPT.RNAData, ::FPTRegressionModel; kwargs...) = param

@testset "First-passage histogram regression" begin
    Q=sparse([-2. 0.;2. 0.]);p=[1.,0.];bins=[.1,.2,.3,.4];edges=[0.;bins]
    truecdf=1 .-exp.(-2edges);truep=diff(truecdf)/last(truecdf)
    @test SGFPT._dwell_bin_edges(bins)==edges
    @test_throws ArgumentError SGFPT._dwell_bin_edges([0.,1.,2.])
    @test_throws ArgumentError SGFPT._dwell_bin_edges([2.,1.])
    @test_throws ArgumentError SGFPT._dwell_bin_edges([1.,Inf])
    @test_throws ArgumentError SGFPT._dwell_bin_edges([1.])
    @test SGFPT._dwell_bin_edges([1.];leftedge=0.)==[0.,1.]
    for method in (1,2),generator in (Q,Matrix(Q))
        values=SGFPT.time_evolve(edges,generator,p,method)
        @test values isa Matrix
        @test size(values)==(5,2)
        @test vec(sum(values;dims=2))≈ones(5) atol=1e-8
        cdf=SGFPT.dwelltimeCDF(bins,generator,[2],p,method)
        mass=SGFPT.dwelltimePDF(bins,generator,[2],p,method)
        @test cdf≈truecdf atol=2e-6
        @test mass≈truep atol=2e-6
        @test length(mass)==length(bins)
    end
    repeated=SGFPT.time_evolve([0.,0.,1.],Q,p,1)
    @test size(repeated)==(3,2)
    @test repeated[1,:]==repeated[2,:]
    @test SGFPT.time_evolve([0.,0.],Q,p,1)==[1. 0.;1. 0.]
    @test_throws ArgumentError SGFPT.dwelltimePDF([0.,1.,2.],Q,[2],p,1)

    # Offset interval + nonexponential law: evolve from actual binding time 0.
    Q3=sparse([-1. 0. 0.;1. -1. 0.;0. 1. 0.]);p3=[1.,0.,0.]
    F(t)=1-exp(-t)*(1+t);later=[2.,3.,4.];lateredges=[1.;later]
    truth=F.(lateredges);mass=diff(truth)/sum(diff(truth))
    for method in (1,2),generator in (Q3,Matrix(Q3))
        @test SGFPT.dwelltimeCDF(later,generator,[3],p3,method)≈truth atol=5e-6
        @test SGFPT.dwelltimePDF(later,generator,[3],p3,method)≈mass atol=5e-6
        irregular=[.7,1.4,3.5];explicit=[.2;irregular]
        expected=diff(F.(explicit));expected./=sum(expected)
        @test SGFPT.dwelltimePDF(irregular,generator,[3],p3,method;leftedge=.2)≈expected atol=5e-6
    end
    @test sum(diff(truecdf))<.56
    @test sum(SGFPT.dwelltimePDF(bins,Q,[2],p,1))≈1.
    # ForwardDiff through the actual ODE/binning path (not just the log guard).
    function first_bin_probability(x)
        z=zero(x[1]);o=one(x[1]);q=sparse([-x[1] z;x[1] z])
        SGFPT.dwelltimePDF(bins,q,[2],[o,z],1)[1]
    end
    grad=ForwardDiff.gradient(first_bin_probability,[2.])[1]
    finite=(first_bin_probability([2.0001])-first_bin_probability([1.9999]))/.0002
    @test grad≈finite rtol=2e-4 atol=1e-6

    # Native TDComponents passes survival rather than absorbed probability.
    transitions=([1,2],[2,1],[2,3],[3,1]);r=[.3,2.,.5,.1,.2,.01]
    types=["ONG","OFFG"];b=[collect(.1:.1:2.),collect(.1:.1:2.)]
    reporter,comp=SGFPT.make_reporter_components_DT(transitions,3,0,0,1,"",[[2,3],[2,3]],types)
    full=SGFPT.make_mat(comp.elementsG,r,3);H=Matrix(full)[[2,3],[2,3]]
    survival(t)=sum(exp(H*t)*[1.,0.])
    truthon=-diff(survival.(vcat(0.,b[1])));truthon./=sum(truthon)
    truthoff=diff(1 .-exp.(-r[1].*vcat(0.,b[2])));truthoff./=sum(truthoff)
    for method in (1,2)
        pred=SGFPT.predictedarray(r,comp,b,reporter,types;method=method)
        @test pred[1]≈truthon atol=2e-5
        @test pred[2]≈truthoff atol=2e-5
        r2=copy(r);r2[1]*=10
        pred2=SGFPT.predictedarray(r2,comp,b,reporter,types;method=method)
        @test pred2[1]≈pred[1] atol=2e-5
    end
    pss=SGFPT.steady_state_vector(full)
    init=SGFPT.init_S(r,[2,3],comp.elementsG,pss,[1,2,3])
    @test sum(init)≈1.
    @test init≈[0.,1.,0.]

    SGFPT.set_onoff([[2]],[[1.,2.,3.]],1)
    @test_throws ArgumentError SGFPT.set_onoff([[2]],[[0.,1.,2.]],1)
    @test_throws ArgumentError SGFPT.set_onoff([[2]],[[2.,3.,4.]],1)
    @test_throws ArgumentError SGFPT.set_onoff([[2]],[[1.,2.,4.]],1)
    hist=zeros(Int,3)
    SGFPT.firstpassagetime!(hist,[0.],[1.],1.5,1.,3,1)
    @test hist==[1,0,0]

    model=FPTRegressionModel()
    data(h)=SGFPT.DwellTimeData("regression","test",[[1.,2.]],[h],["ONG"])
    for ll in (SGFPT.loglikelihood,SGFPT.loglikelihood_ad)
        @test first(ll([.2,.8],data([2.,3.]),model))≈2log(.2)+3log(.8)
        @test first(ll([-.2,1.2],data([1.,1.]),model))==-Inf
        @test first(ll([0.,1.],data([1.,3.]),model))==-Inf
        @test first(ll([0.,1.],data([0.,4.]),model))==0.
        @test first(ll([1e-30,1.],data([1.,3.]),model))≈log(1e-30)
        @test first(ll([NaN,1.],data([1.,3.]),model))==-Inf
    end
    gradient=ForwardDiff.gradient(x->first(SGFPT.loglikelihood_ad(x,data([2.,3.]),model)),[.2,.8])
    @test gradient≈[2/.2,3/.8]
    @test !SGFPT._histogram_predictions_invalid([-1e-15,1.],[0.,1.])
    @test !SGFPT._histogram_predictions_invalid([-4e-12,1.],[0.,1.])
    @test SGFPT._histogram_predictions_invalid([-4e-12,1.],[1.,1.])
    @test SGFPT._histogram_predictions_invalid([-1e-5,1.],[0.,1.])

    # The guard is shared with RNA histograms, not just first passage data.
    rna=SGFPT.RNAData("regression","test",2,[2.,3.],1.)
    for ll in (SGFPT.loglikelihood,SGFPT.loglikelihood_ad)
        @test first(ll([.2,.8],rna,model))≈2log(.2)+3log(.8)
        @test first(ll([0.,1.],rna,model))==-Inf
        @test first(ll([1e-30,1.],rna,model))≈2log(1e-30)
    end
    @test only(Zygote.gradient(x->first(SGFPT.loglikelihood_ad(x,rna,model)),[.2,.8]))≈[10.,3/.8]
end
