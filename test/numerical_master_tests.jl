using Test, LinearAlgebra, SparseArrays, Random
using StochasticGene

@testset "Numerical master-equation extensions" begin
    SG = StochasticGene
    tr = ([1,2], [2,1])
    comp = SG.MComponents(tr, 2, 0, 45, 1.0, "")
    rates = [.3, .5, 8.]
    prob = transient_master_problem(comp, rates)
    P0 = transient_master_initial([.625,.375], 45)
    reference = reshape(exp(Matrix(SG.make_mat_M(comp,rates))) * vec(P0), size(P0))

    @testset "Independent full-CME references" begin
        rk = transient_master_closure(prob,P0,1.;ode_steps=200)
        tay = transient_master_closure_taylor(prob,P0,1.;n_steps=40)
        et = transient_master_closure_exp_taylor(prob,P0,1.;n_steps=20)
        for actual in (rk,tay,et)
            @test actual ≈ reference atol=2e-9
            @test sum(actual) ≈ 1 atol=1e-10
            @test minimum(actual) >= -1e-12
        end
        @test transient_master_closure(prob,P0,1.;return_joint=false) ≈ vec(sum(rk;dims=1))
        rich = transient_master_strang_richardson(prob,P0,1.,25)
        coarse = transient_master_strang(prob,P0,1.,25)
        @test norm(rich-reference,1) < norm(coarse-reference,1)/100
        @test transient_master_strang_purebd(prob,P0,1.,200) ≈ reference atol=1e-5

        rc = SG.MComponents(tr,2,1,45,1.,"")
        rr = [.3,.5,2.,8.]
        general = transient_master_problem_general(rc,rr)
        init = transient_master_initial(fill(1/rc.nT,rc.nT),45)
        full = reshape(exp(Matrix(SG.make_mat_M(rc,rr))) * vec(init),size(init))
        @test_throws ArgumentError transient_master_problem(rc,rr)
        for solve in (transient_master_closure, transient_master_closure_taylor, transient_master_closure_exp_taylor)
            @test solve(general,init,1.) ≈ full atol=1e-8
        end
        @test transient_master_strang_purebd(general,init,1.,200) ≈ full atol=1e-5
    end

    @testset "Window semantics, zero decay, and stable thinning" begin
        poisson = SG.TransientMasterProblem(zeros(1,1),[4.],0.)
        initial = transient_master_initial([1.],31)
        p = transient_master_closure(poisson,initial,1.)
        expected = [exp(-4.)*4.0^m/factorial(big(m)) for m in 0:30]
        @test vec(p) ≈ Float64.(expected) atol=1e-9
        @test transient_master_strang(poisson,initial,1.,10) ≈ p atol=1e-9
        @test transient_master_closure(poisson,initial,1.;M_out=3) ≈ p[:,1:4] atol=1e-10
        @test sum(transient_master_closure(poisson,initial,1.;M_out=3)) < .5
        @test transient_master_closure(prob,P0,0.;M_out=4) == P0[:,1:5]
        high = zeros(1,1001);high[1,end]=1.
        thin = SG._compose_initial_after_death(high,.99,1000)
        @test sum(thin) ≈ 1 atol=2e-11
        @test sum((0:1000).*vec(thin)) ≈ 990 atol=1e-8
        small = SG._compose_initial_after_death(high,.1,110)
        @test small ≈ SG._compose_initial_after_death(high,.1,1000)[:,1:111]
        pois = SG._poisson_convolve([1.;zeros(1200)],1000.)
        @test sum(pois) ≈ 1 atol=1e-8
        @test argmax(pois) in (1000,1001)
    end

    @testset "Stiff phi functions" begin
        L = [-1000. 1000.;1000. -1000.]
        phi = SG._compute_phi_matrices(L,1.,8)
        @test phi[1] ≈ exp(L) atol=1e-12
        @test phi[1] ≈ fill(.5,2,2) atol=1e-12
        for j in 1:8
            @test L*phi[j+1] ≈ phi[j]-Matrix{Float64}(I,2,2)/factorial(j-1) atol=2e-11
        end
        stiff = SG.TransientMasterProblem(L,[.1,.2],1.)
        initial = transient_master_initial([1.,0.],12)
        @test_throws ArgumentError transient_master_closure_exp_taylor(stiff,initial,1.;n_steps=10)
        actual = transient_master_closure_exp_taylor(stiff,initial,1.;n_steps=2000)
        M = SG.make_mat_M(sparse(L),spdiagm(0=>[.1,.2]),1.,12)
        full = reshape(exp(Matrix(M))*vec(initial),size(initial))
        @test all(isfinite,actual)
        @test actual ≈ full atol=1e-7
    end

    @testset "State-dependent decay" begin
        A,B = SG._problem_AB(prob)
        mus = [0.,1.5]
        stateprob = TransientMasterProblemStateMu(A,B,mus)
        nT,nh=size(P0)
        M=kron(Matrix{Float64}(I,nh,nh),A)
        for m in 0:nh-1
            rows=m*nT+1:(m+1)*nT
            M[rows,rows] .-= m.*Diagonal(mus)
            m>0 && (M[(m-1)*nT+1:m*nT,rows] .+= m.*Diagonal(mus))
            m<nh-1 && (M[(m+1)*nT+1:(m+2)*nT,rows] .+= B)
        end
        full=reshape(exp(M)*vec(P0),size(P0))
        @test transient_master_strang_statewise(stateprob,P0,1.,200) ≈ full atol=1e-5
        shifted80=transient_master_strang_shifted(stateprob,P0,1.,80)
        shifted160=transient_master_strang_shifted(stateprob,P0,1.,160)
        @test shifted160 ≈ full atol=1e-5
        @test norm(shifted160-full) < norm(shifted80-full)/3.5
        @test_throws ArgumentError transient_master_strang_shifted(stateprob,P0,1.,10;mu0=2.)
        rc=SG.MComponents(tr,2,1,20,1.,"")
        gen=transient_master_problem_general(rc,[.3,.5,2.,8.])
        sp=TransientMasterProblemStateMu(gen.A,gen.B,[0.,.5,1.,1.5])
        init=transient_master_initial(fill(1/rc.nT,rc.nT),20)
        @test_throws ArgumentError transient_master_strang_statewise(sp,init,1.,10)
        @test all(isfinite,transient_master_strang_shifted(sp,init,.1,10))
    end

    @testset "Block-Thomas matches current capped generators" begin
        rng=MersenneTwister(717)
        for R in (0,1,2), nhist in (2,12,30), _ in 1:2
            c=SG.MComponents(tr,2,R,nhist,1.,"")
            r=exp.(randn(rng, R==0 ? 3 : R+3))
            M=SG.make_mat_M(c,r)
            p=SG.steady_state_vector(M;solver=:closure_thomas,nT=c.nT)
            ref=SG.steady_state_vector(M;solver=:augmented)
            @test p ≈ ref atol=1e-9
            @test sum(p) ≈ 1 atol=1e-12
            @test minimum(p)>=0
            @test norm(M*p,Inf)<1e-9
            @test SG.steady_state(M,c.nT,2;steady_state_solver=:closure_thomas) ≈ SG.steady_state(M,c.nT,2) atol=1e-9
        end
        M=SG.make_mat_M(comp,rates)
        @test_throws ArgumentError SG.steady_state_vector(M;solver=:closure_thomas)
        @test_throws ArgumentError SG.steady_state_vector(M;solver=:closure_thomas,nT=0)
        @test_throws ArgumentError SG.steady_state_vector(M;solver=:closure_thomas,nT=7)
        bad=copy(M);bad[1,5]+=.1
        @test_throws ArgumentError SG.steady_state_vector(bad;solver=:closure_thomas,nT=2)
        burst=SG.MComponents(tr,2,0,20,1.,"",2)
        @test_throws ArgumentError transient_master_problem(burst,rates)
        @test_throws ArgumentError transient_master_problem_general(burst,rates)
        @test_throws ArgumentError SG.steady_state_vector(SG.make_mat_M(burst,rates);solver=:closure_thomas,nT=2)
        no_birth=SG.make_mat_M(comp,[.3,.5,0.])
        pzero=SG.steady_state_vector(no_birth;solver=:closure_thomas,nT=2)
        @test sum(pzero[1:2]) ≈ 1 atol=1e-10
        @test norm(no_birth*pzero,Inf)<1e-10
        cap1=sparse([-.3 .5;.3 -.5])
        @test SG.steady_state_vector(cap1;solver=:closure_thomas,nT=2) ≈ [.625,.375]
        @test isfinite(SG.ForwardDiff.derivative(8.) do x
            SG.steady_state_vector(SG.make_mat_M(comp,[.3,.5,x]);solver=:closure_thomas,nT=2)[1]
        end)
    end

    @testset "Input checks" begin
        @test_throws ArgumentError SG.TransientMasterProblem(zeros(1,1),[1.],-1.)
        @test_throws ArgumentError TransientMasterProblemAB(zeros(1,1),ones(1,1),1.)
        @test_throws ArgumentError transient_master_closure(prob,P0,-1.)
        @test_throws ArgumentError transient_master_closure(prob,P0,0.;ode_steps=0)
        @test_throws ArgumentError transient_master_closure_taylor(prob,P0,1.;R=0)
        @test_throws ArgumentError transient_master_closure_exp_taylor(prob,P0,1.;R=21)
        @test_throws ArgumentError transient_master_strang(prob,P0,[-1.,0.])
        @test_throws DimensionMismatch transient_master_closure(prob,ones(1,4),1.)
    end
end
