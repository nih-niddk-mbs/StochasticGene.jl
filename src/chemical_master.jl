"""
offonPDF(t::Vector,r::Vector,n::Int,nr::Int)

Active (ON) and Inactive (OFF) time distributions for GRSM model
Takes difference of ON and OFF time CDF to produce PDF
"""
function offonPDF(t::Vector,r::Vector,n::Int,nr::Int,method)
    T,TA,TI = mat_GSR_DT(r,n,nr)
    pss = normalized_nullspace(T)
    SA=ontimeCDF(t,r,n,nr,TA,pss,method)
    SI=offtimeCDF(t,r,n,nr,TI,pss,method)
    return pdf_from_cdf(t,SI), pdf_from_cdf(t,SA)
end
function offonPDF_offeject(t::Vector,r::Vector,n::Int,nr::Int,method)
    T,TA,TI = mat_GSR_DT_offeject(r,n,nr)
    pss = normalized_nullspace(T)
    SA=ontimeCDF(t,r,n,nr,TA,pss,method)
    SI=offtimeCDF(t,r,n,nr,TI,pss,method)
    return pdf_from_cdf(t,SI), pdf_from_cdf(t,SA)
end
# function offPDF(t::Vector,r::Vector,n::Int,nr::Int)
#     T,_,TI = mat_GSR(r,n,nr)
#     pss = normalized_nullspace(T)
#     SI=offtimeCDF(t,r,n,nr,TI,pss)
#     pdf_from_cdf(t,SI)
# end

"""
pdf_from_cdf(t,S)
"""
function pdf_from_cdf(t,S)
    P = diff(S)
    P/(sum(P)*(t[2]-t[1]))
end
"""
onCDF(t::Vector,r::Vector,n::Int,nr::Int)

"""
# function onCDF(t::Vector,r::Vector,n::Int,nr::Int)
#     T,TA,TI = mat_GSR(r,n,nr)
#     pss = normalized_nullspace(T)
#     ontimeCDF(t,r,n,nr,TA,pss)
# end
"""
ontimeCDF(tin::Vector,n::Int,nr::Int,rin::Vector,TA::Matrix,pss::Vector)
offtimeCDF(tin::Vector,n::Int,nr::Int,r::Vector,TI::Matrix,pss::Vector)

ON(OFF) dwell time distributions of GRS model
Found by computing accumulated probability into OFF(ON) states
where transitions out of OFF(ON) states are zeroed, starting from first instance of ON(OFF) state
weighted by steady state distribution (i.e. solving a first passage time problem)x
"""
function ontimeCDF(tin::Vector,rin::Vector,n::Int,nr::Int,TA::Matrix,pss::Vector,method)
    t = [tin ; tin[end] + tin[2]-tin[1]] #add a time point so that diff() gives original length
    SAinit = init_prob(pss,n,nr)
    SA = time_evolve(t,TA,SAinit,method)
    accumulate(SA,n,nr)  # accumulated prob into OFF states
end
function offtimeCDF(tin::Vector,r::Vector,n::Int,nr::Int,TI::Matrix,pss::Vector,method)
    t = [tin ; tin[end] + tin[2]-tin[1]]
    nonzerosI = nonzero_rows(TI)  # only keep nonzero rows to reduce singularity of matrix
    TI = TI[nonzerosI,nonzerosI]
    SIinit = init_prob(pss,r,n,nr,nonzerosI)
    SI = time_evolve(t,TI,SIinit,method)
    accumulate(SI,n,nr,nonzerosI,length(pss)) # accumulated prob into ON states
end
"""
steady_state_offpath(rin::Vector,n::Int,nr::Int,nhist::Int,nalleles::Int)
GRS model where mRNA decay rate is accelerated to account for nonviability of off-pathway pre-mRNA
from RNA that is recursively spliced
"""
function steady_state_offpath(rin::Vector,n::Int,nr::Int,nhist::Int,nalleles::Int)
    r = copy(rin)
    nu = get_nu(r,n,nr)
    eta = get_eta(r,n,nr)
    # r[end] /= survival_fraction(nu,eta,nr)
    yieldfactor = survival_fraction(nu,eta,nr)
    mhist=steady_state(r,n,nr,nhist,nalleles)
    technical_loss(mhist,yieldfactor,nhist)
end
"""
steady_state(rin::Vector,n::Int,nr::Int,nhist::Int,nalleles::Int)
Steady state distribution of mRNA in GRM model (which is the same as GRSM model)
"""
function steady_state(rin::Vector,n::Int,nr::Int,nhist::Int,nalleles::Int)
    r = rin/rin[end]
    gammap,gamman = get_gamma(r,n)
    nu = get_nu(r,n,nr)
    T,B = transition_rate_mat(n,nr,gammap,gamman,nu)
    P = initial_pmf(T,nu[end],n,nr,nhist)
    mhist=steady_state(nhist,P,T,B)
    allele_convolve(mhist[1:nhist],nalleles) # Convolve to obtain result for n alleles
end
"""
steady_state_offeject(rin::Vector,n::Int,nr::Int,nhist::Int,nalleles::Int)
Steady state distribution of mRNA in GRM model (which is the same as GRSM model)
"""
function steady_state_offeject(rin::Vector,n::Int,nr::Int,nhist::Int,nalleles::Int)
    r = rin/rin[end]
    gammap,gamman = get_gamma(r,n)
    nu = get_nu(r,n,nr)
    eta = get_eta(r,n,nr)
    T,B = transition_rate_mat_offeject(n,nr,gammap,gamman,nu,eta)
    P = initial_pmf(T,nu[end],n,nr,nhist)
    mhist=steady_state(nhist,P,T,B)
    allele_convolve(mhist[1:nhist],nalleles) # Convolve to obtain result for n alleles
end
"""
steady_state(nhist::Int,nalleles::Int,ejectrate,P,T,B,tol = 1e-6)
Iterative algorithm for computing null space of truncated transition rate matrix
of Master equation of GR model to give steady state of mRNA in GRM model
for single allele
"""
function steady_state(nhist::Int,P::Matrix,T::Matrix,B::Matrix,tol = 1e-6,stepmax=1000)
    total = size(P,2)
    steps = 0
    err = 1.
    A = T - B
    while err > tol && steps < stepmax
        P0 = copy(P)
        P[:,1] = try -A\P[:,2]
        catch
            P[:,1] = (-A + UniformScaling(1e-18))\P[:,2]
        end
        for m = 2:total-1
            P[:,m] = @inbounds -(A - UniformScaling(m-1))\((B*P[:,m-1]) + m*P[:,m+1])
        end
        P[:,total] = -(T - UniformScaling(total-1))\((B*P[:,total-1]))
        P /=sum(P)
        err = norm(P-P0,Inf)
        steps += 1
    end
    sum(P,dims=1)   # marginalize over GR states
end
"""
steady_state(r,n,nhist,nalleles)
Steady State of mRNA in G (telelgraph) model
"""
function steady_state(r::Vector,yieldfactor::Float64,n::Int,nhist::Int,nalleles::Int)
    mhist = steady_state(r,n,nhist_loss(nhist,yieldfactor),nalleles)
    technical_loss(mhist,yieldfactor,nhist)
end

function steady_state(r::Vector,n::Int,nhist::Int,nalleles::Int)
    P = steady_state_full(r,n,nhist)
    steady_state_rna(P,n,nhist,nalleles)
end

function steady_state_full(r::Vector,n::Int,nhist::Int)
    M = mat_GM(r,n,nhist)
    normalized_nullspace(M)
end

function steady_state_rna(P,n,nhist,nalleles)
    mhist = marginalize(P,n)
    allele_convolve(mhist,nalleles)[1:nhist]
end
"""
nhist_loss(nhist,yieldfactor)
"""
nhist_loss(nhist,yieldfactor) = round(Int,nhist/yieldfactor)

"""
mat(r,n,nr,nhist)
Transition rate matrices for GSR model
return T,TA,TI
"""
function mat_GSR_DT(r,n,nr)
    gammap,gamman = get_gamma(r,n)
    nu = get_nu(r,n,nr)
    eta = get_eta(r,n,nr)
    T = transition_rate_mat_T(n,nr,gammap,gamman,nu,eta)
    transition_rate_mat(T,n,nr)
end
function mat_GSR_DT_offeject(r,n,nr)
    gammap,gamman = get_gamma(r,n)
    nu = get_nu(r,n,nr)
    eta = get_eta(r,n,nr)
    T,_ = transition_rate_mat_offeject(n,nr,gammap,gamman,nu,eta)
    transition_rate_mat(T,n,nr,2)
end

function mat_GSR_T(r,n,nr)
    gammap,gamman = get_gamma(r,n)
    nu = get_nu(r,n,nr)
    eta = get_eta(r,n,nr)
    transition_rate_mat_T(n,nr,gammap,gamman,nu,eta)
end

function get_rates_GSR(r,n,nr)
    gammap,gamman = get_gamma(r,n)
    nu = get_nu(r,n,nr)
    eta = get_eta(r,n,nr)
    return gammap,gamman,nu,eta
end

"""
mat_GM(r,n,nhist)
Transition rate matrix of GM model
"""
function mat_GM(r,n,nhist)
    gammap,gamman = get_gamma(r,n)
    transition_rate_mat(n,gammap,gamman, r[2*n+1],r[2*n+2],nhist)
end
"""
transient(ts::Vector,r,n,nhist,nalleles,P0)

Compute mRNA pmf of GM model at times in vector ts starting
with initial condition P0
"""

"""
transient(t,r,yieldfactor,n,nhist,nalleles,P0::Vector,method)
transient(t,n,nhist,nalleles,P0,Mvals,Mvects)

Compute mRNA pmf of GM model at time t given initial condition P0
and eigenvalues and eigenvectors of model transition rate matrix
method = 1 uses JuliaDifferentialEquations.jl
method != 1 uses eigenvalue decomposition
"""
function transient(t,r,yieldfactor,n,nalleles,P0::Vector,method)
    mhist = transient(t,r,n,nalleles,P0,method)
    technical_loss(mhist,yieldfactor)
end
function transient(t::Float64,r,n,nalleles,P0::Vector,method)
    P = transient(t,r,n,P0,method)
    mh = marginalize(P,n)
    allele_convolve(mh,nalleles)
end
function transient(t::Vector,r,n,nalleles,P0::Vector,method)
    P = transient(t,r,n,P0,method)
    mhist = Array{Array,1}(undef,length(t))
    for i in 1:size(P,1)
        mh = marginalize(P[i,:],n)
        mhist[i]= allele_convolve(mh,nalleles)
    end
    return mhist
end
function transient(t,r,n,P0::Vector,method)
    nhist = Int(length(P0)/(n+1)) - 2
    M = mat_GM(r,n,nhist)
    time_evolve(t,M,P0,method)
end

function transient_delay(t::Vector,r0::Vector,r1::Vector,delay::Float64,n::Int,nalleles,P0::Vector)
    P = time_evolve_delay(t,r0,r1,delay,n,P0)
    mhist = Array{Array,1}(undef,length(t))
    for i in 1:size(P,1)
        mh = marginalize(P[i,1:end-1],n)
        mhist[i]= allele_convolve(mh,nalleles)
    end
    return mhist
end


"""
time_evolve(t,M::Matrix,Sinit::Vector)
Eigenvalue solution of Linear ODE with rate matrix T and initial vector Sinit
"""
function time_evolve(t,M::Matrix,S0::Vector,method)
    if method == 1
        return time_evolve_diff(t,M,S0)
    else
        return time_evolve(t,M,S0)
    end
end

function time_evolve(t,M::Matrix,Sinit::Vector)
    vals,vects = eig_decompose(M)
    weights = solve_vector(vects,Sinit)
    time_evolve(t,vals,vects,weights)
end
"""
time_evolve(t,vals::Vector,vects::Matrix,weights::Vector)

"""
function time_evolve(t::Float64,vals::Vector,vects::Matrix,weights::Vector)
    n = length(vals)
    S = zeros(n)
    for j = 1:n
        for i = 1:n
            S[j] += real(weights[i]*vects[j,i]*exp.(vals[i]*t))
        end
    end
    return S
end
function time_evolve(t::Vector,vals::Vector,vects::Matrix,weights::Vector)
    ntime = length(t)
    n = length(vals)
    S = Array{Float64,2}(undef,ntime,n)
    for j = 1:n
        Sj = zeros(ntime)
        for i = 1:n
            Sj += real(weights[i]*vects[j,i]*exp.(vals[i]*t))
        end
        S[:,j]=Sj
    end
    return S
end
"""
Solve transient problem using DifferentialEquations.jl
"""
function time_evolve_diff(t,M::Matrix,P0)
    global M_global = copy(M)
    tspan = (0.,t[end])
    prob = ODEProblem(fglobal,P0,tspan)
    # sol = solve(prob,saveat=t, lsoda(),abstol = 1e-4, reltol = 1e-4)
    sol = solve(prob,saveat=t, lsoda())
    return sol'
end

fglobal(u,p,t) = M_global*u

function time_evolve_delay(t,r0,r1,delay,n,P0)
    tspan = (0.,t[end])
    nhist = Int(length(P0)/(n+1)) - 2
    p  = [r0;r1;delay;n;nhist]
    P0 = [P0;1.]
    prob = ODEProblem(fdelay!,P0,tspan,p)
    sol = solve(prob,saveat=t)
    return sol'
end

function fdelay!(du,u,p,t)
    n = Int(p[end-1])
    nhist = Int(p[end])
    delay = p[end-3]
    r0 = p[1:2*(n+1)]
    r1 = p[2*(n+1)+1:4*(n+1)]
    r = r0*u[end] + (1-u[end])*r1
    M = mat_GM(r,n,nhist)
    du[1:end-1] = M*u[1:end-1]
    du[end] = -delay*u[end]
end


"""
technical_loss(mhist,yieldfactor,nhist)

"""
function technical_loss!(mhist,yieldfactor)
    for i in eachindex(mhist)
        mhist[i] = technical_loss(mhist[i],yieldfactor,length(mhist[i]))
    end
end
function technical_loss_poisson(mhist,yieldfactor,nhist)
    p = zeros(nhist)
    for m in eachindex(mhist)
        d = Poisson(yieldfactor*(m-1))
        for c in 1:nhist
            p[c] += mhist[m]*pdf(d,c-1)
        end
    end
    normalize_histogram(p)
end
function technical_loss(mhist,yieldfactor,nhist)
    p = zeros(nhist)
    for m in eachindex(mhist)
        d = Binomial(m-1,clamp(yieldfactor,0.,1.))
        for c in 1:nhist
            p[c] += mhist[m]*pdf(d,c-1)
        end
    end
    normalize_histogram(p)
end

function additive_noise(mhist,noise,nhist)
    p = zeros(nhist)
    d = Poisson(noise)
    for m in 1:nhist
        for n in 1:m
            p[m] += mhist[n]*pdf(d,m-n)
        end
    end
    normalize_histogram(p)
end


function threshold_noise(mhist,noise,yieldfactor,nhist)
    h = additive_noise(mhist,noise,nhist)
    technical_loss(h,yieldfactor,nhist)
end

"""
solve_vector(A::Matrix,b::vector)
solve A x = b
If matrix divide has error higher than tol
use SVD and pseudoinverse with threshold
"""
function solve_vector(A::Matrix,b::Vector,th = 1e-16,tol=1e-1)
    x = A\b
    if norm(b-A*x,Inf) > tol
        M = svd(A)
        Sv = M.S
        Sv[abs.(Sv) .< th] .= 0.
        Sv[abs.(Sv) .>= th] = 1 ./ Sv[abs.(Sv) .>= th]
        x = M.V * diagm(Sv) * M.U' * b
    end
    return x[:,1] # return as vector
end
"""
initial_pmf(T,ejectrate,n,nr,nhist)

Tensor product of nullspace of T and Poisson density
with rate = mean ejection rate of mRNA
Used for initial condition of iteration algorithm to
compute steady state pmf of GRM Master Equation
"""
function initial_pmf(T,ejectrate,n,nr,nhist)
    t0 = normalized_nullspace(T)
    ejectrate *= sum(t0[(n+1)*2^(nr-1)+1:end])
    d = Poisson(ejectrate)
    nT = length(t0)
    total = nhist + 2
    P = Array{Float64,2}(undef,nT,total)
    for m = 1:total
        P[:,m] = t0 * pdf(d,m)
    end
    P
end
"""
get_gamma(r,n,nr)
G state forward and backward transition rates
for use in transition rate matrices of Master equation
(different from gamma used in Gillespie algorithms)
"""
function get_gamma(r,n)
    gammaf = zeros(n+2)
    gammab = zeros(n+2)
    for i = 1:n
        gammaf[i+1] = r[2*(i-1)+1]
        gammab[i+1] = r[2*i]
    end
    return gammaf, gammab
end
"""
get_nu(r,n,nr)
R step forward transition rates
"""
function get_nu(r,n,nr)
    r[2*n+1 : 2*n+nr+1]
end
"""
get_eta(r,n,nr)
Intron ejection rates at each R step
"""
function get_eta(r,n,nr)
    eta = zeros(nr)
    if length(r) > 2*n + 2*nr
        eta[1] = r[2*n + 1 + nr + 1]
        for i = 2:nr
            eta[i] = eta[i-1] + r[2*n + 1 + nr + i]
        end
    end
    return eta
end
"""
survival_fraction(nu,eta,nr)
Fraction of introns that are not spliced prior to ejection
"""
function survival_fraction(nu,eta,nr)
    pd = 1.
    for i = 1:nr
        pd *= nu[i+1]/(nu[i+1]+eta[i])
    end
    return pd
end

function reparametrize(r,n)
    gammaf,gammb = get_gamma(r,n)
    alpha = Vector{Float64}(undef,n)
    alpha[1] = 1
    for i in 2:n
        alpha[i] = alpha[i-1] * gammaf[i]/gammab[i]
    end
    return alpha
end

function invert_parameter(alpha,n)


end

"""
normalized_nullspace(M::AbstractMatrix)
Compute the normalized null space of a nxn matrix
of rank n-1 using QR decomposition with pivoting
"""
function normalized_nullspace(M::AbstractMatrix)
    F = qr(M,Val(true));  #QR decomposition with pivoting
    R = F.R
    m = size(M,1)
    # if rank(M) == m-1
    p = zeros(m)
    # Back substitution to solve R*p = 0
    p[end] = 1.
    for i in 1:m-1
        p[m-i] = -R[m-i,m-i+1:end]'*p[m-i+1:end]/R[m-i,m-i]
    end
    p /= sum(p);
    return F.P*p
    # else
    #     return zeros(m)
    # end
end
"""
allele_convolve(mhist,nalleles)
Convolve to compute distribution for contributions from multiple alleles
"""
function allele_convolve(mhist,nalleles)
    nhist = length(mhist)
    mhists = Array{Array{Float64,1}}(undef,nalleles)
    mhists[1] = float.(mhist)
    for i = 2:nalleles
        mhists[i] = zeros(nhist)
        for m = 0:nhist-1
            for m2 = 0:min(nhist-1,m)
                mhists[i][m+1] += mhists[i-1][m-m2+1]*mhist[m2+1]
            end
        end
    end
    return mhists[nalleles]
end
"""
init_prob(pss,n,nr)
Initial condition for first passage time calculation of active (ON) state
(Sum over all states with first R step occupied weighted by equilibrium state
of all states with all R steps not occupied with reporters)
"""
function init_prob(pss,n,nr)
    l = length(pss)
    base = findbase(l,n,nr)
    SAinit = zeros(l)
    if base == 3
        for z=1:2^(nr-1)
            aSAinit = (n+1) + (n+1)*decimal(vcat(2,digits(z-1,base=2,pad=nr-1)))
            apss = (n+1) + (n+1)*decimal(vcat(0,digits(z-1,base=2,pad=nr-1)))
            SAinit[aSAinit] = pss[apss]
        end
    else
        for z=1:2^(nr-1)
            aSAinit = (n+1) + (n+1)*decimal(vcat(1,zeros(Int,nr-1)))
            apss = n+1
            SAinit[aSAinit] = pss[apss]
        end

    end
    SAinit / sum(SAinit)
end
"""
init_prob(pss,n,nr)
Initial condition for first passage time calculation of inactive (OFF) state
(Sum over all states with no R step introns occupied weighted by equilibrium state
of all states with a single R step intron not)
"""
function init_prob(pss,r,n,nr,nonzeros)
    l = length(pss)
    base = findbase(l,n,nr)
    SIinit = zeros(l)
    nu = get_nu(r,n,nr)
    eta = get_eta(r,n,nr)
    if base == 3
        # Start of OFF state by ejection
        for i in 1:n+1, z in 1:2^(nr-1)
            exons = digits(z-1,base=2,pad=nr-1)
            ainit = i + (n+1)*decimal(vcat(exons,0))
            apss = i + (n+1)*decimal(vcat(exons,2))
            SIinit[ainit] += pss[apss]*nu[nr+1]
        end
        # Start of OFF state by splicing
        for i in 1:n+1, z in 1:2^(nr)
            exons = digits(z-1,base=2,pad=nr)
            intronindex = findall(exons.==1)
            ainit = i + (n+1)*decimal(exons)
            for j in intronindex
                introns = copy(exons)
                introns[j] = 2
                apss  = i + (n+1)*decimal(introns)
                SIinit[ainit] += pss[apss]*eta[j]
            end
        end
    else
        # Start of OFF state by ejection
        for i in 1:n+1
            ainit = i
            apss = i + (n+1)*decimal(vcat(zeros(Int,nr-1),1),2)
            SIinit[ainit] += pss[apss]*nu[nr+1]
        end
        # Start of OFF state by splicing
        for i in 1:n+1
            ainit = i
            for j in 1:nr
                introns = zeros(Int,nr)
                introns[j] = 1
                apss  = i + (n+1)*decimal(introns,2)
                SIinit[ainit] += pss[apss]*eta[j]
            end
        end

    end
    SIinit = SIinit[nonzeros]
    SIinit/sum(SIinit)
end
"""
accumulate(SA::Matrix,n,nr)
Sum over all probability vectors accumulated into OFF states
"""
function accumulate(SA,n,nr)
    l,p = size(SA)
    base = findbase(p,n,nr)
    SAj = zeros(size(SA)[1])
    for i=1:n+1, z=1:base^nr
        zdigits = digits(z-1,base=base,pad=nr)
        if ~any(zdigits.> base-2)
            a = i + (n+1)*(z-1)
            SAj += SA[:,a]
        end
    end
    return SAj
end
"""
accumulate(SI::Matrix,n,nr,nonzeros)
Sum over all probability vectors accumulated into ON states
"""
function accumulate(SI,n,nr,nonzeros,p)
    # Sum over all probability vectors accumulated into ON states
    l = size(SI)[1]
    base = findbase(p,n,nr)
    SIj = zeros(size(SI)[1])
    for i=1:n+1, z=1:base^nr
        zdigits = digits(z-1,base=base,pad=nr)
        if any(zdigits.> base-2)
            a = i + (n+1)*(z-1)
            if a in nonzeros
                SIj += SI[:,findfirst(a .== nonzeros)]
            end
        end
    end
    return SIj
end
"""
marginalize(p,n,nhist)
Marginalize over G states
"""
function marginalize(p,n,nhist)
    mhist = zeros(nhist)
    nT = n+1
    for m in 1:nhist
        i = (m-1)*nT
        mhist[m] = sum(p[i+1:i+nT])
    end
    return mhist
end

function marginalize(p,n)
    nhist = Int(length(p)/(n+1)) - 2
    marginalize(p,n,nhist)
end
