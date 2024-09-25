# This file is part of StochasticGene.jl    

# utilities.jl


"""
    invert_dict(D)

invert a dictionary
"""
invert_dict(D) = Dict(D[k] => k for k in keys(D))


"""
    eig_decompose(M)

Take matrix M and return values and vectors
"""
function eig_decompose(M)
    Meiv = eigen(M)
    return Meiv.values, Meiv.vectors
end

"""
    nonzero_rows(T)

Returns an array of row indices that have at least one nonzero element for matrix T
"""
function nonzero_rows(T)
    n = size(T)[2]
    nonzeros = Array{Int}(undef, 0)
    for a in 1:size(T)[1]
        if T[a, :] != zeros(n)
            push!(nonzeros, a)
        end
    end
    return nonzeros
end

"""
    nonzero_states(states,nonzeros)

return reindexed state vector with zeros removed
"""
nonzero_states(states, nonzeros) = [findfirst(o .== nonzeros) for o in intersect(states, nonzeros)]

"""
normalize_histogram(hist)

Returns normalized histogram hist
"""
normalize_histogram(hist) = hist / sum(hist)


"""
var_update(vartuple, newValue)

vartuple = count, mean, M2
iteration step for online algorithm to compute
mean and variance
returns count, mean, M2
"""
function var_update(vartuple, newValue)
    count, mean, M2 = vartuple
    count += 1
    delta = newValue - mean
    mean .+= delta / count
    delta2 = newValue - mean
    M2 .+= delta .* delta2
    return count, mean, M2
end
"""
cov_update(covtuple,newValue1,newValue2)

update covtuple for online algorithm to compute
    covariance
"""
function cov_update(covtuple, newValue1, newValue2)
    count, meanx, meany, C = covtuple
    count += 1
    dx = newx - meanx
    meanx .+= dx / count
    meany .+= (newy - meany) / count
    C .+= dx * (newy - meany)
    return count, meanx, meany, C
end
"""
finalize(vartuple)

Retrieve the mean, variance and sample variance from an aggregate
"""
function finalize(vartuple)
    count, mean, M2 = vartuple
    if count < 2
        return mean, M2
    else
        return mean, M2 ./ (count - 1)
    end
end

"""
    truncate_histogram(x::Array, yield::Float64, nhistmax::Int)

return histogram x that accounts for fraction yield of total sum
"""
function truncate_histogram(x::Array, yield::Float64, nhistmax::Int)
    if yield < 1.0
        counts = sum(x)
        if counts > 0
            nhist = 1
            ratio = 0.0
            # Only keep yield of all data to limit support of histogram
            while ratio <= yield && nhist <= nhistmax
                ratio = sum(x[1:nhist]) / counts
                nhist += 1
            end
            return x[1:min(nhist, length(x))]
        else
            return 0
        end
    else
        return x[1:min(nhistmax, length(x))]
    end
end

"""
combine_histogram(x1::Array,x2::Array)

add histograms x1 and x2

"""
function combine_histogram(x1::Array, x2::Array)
    n = min(length(x1), length(x2))
    x1[1:n] + x2[1:n]
end

"""
pooled_mean(means,counts)

compute weighted average of means given count totals
"""
pooled_mean(means, counts) = counts' * means / sum(counts)

"""
pooled_variance(vars,counts)

computed weighted average of variances given counts
"""
pooled_variance(vars, counts) = (counts .- 1)' * vars / sum(counts .- 1)

"""
pooled_std(std,counts)

compute weighted average of standard deviations given counts
"""
pooled_std(std, counts) = sqrt(pooled_variance(std .^ 2, counts))

"""
var_ratio(mua,mub,vara,varb,cov)

Compute ratio of variances of two variables given means, variances, and covariance of variables
"""
var_ratio(mua, mub, vara, varb, cov) = mua^2 / mub^2 * (vara / mua^2 - 2 * cov / (mua * mub) + varb / mub^2)

"""
decimal(x::Vector)

convert number to base 10 from any base
"""
function decimal(x::Vector, base::Int=3)
    nr = length(x)
    x' * base .^ collect(0:nr-1)
end

"""
    digit_vector(z,base,R)

convert z to vector of length R in base=base
"""
digit_vector(z, base, R) = digits(z - 1, base=base, pad=R)

"""
findbase(l,n,nr)

Find the number of G states for transition rate matrix of size l
"""
function findbase(l, n, nr)
    if l == (n + 1) * 3^nr
        return 3
    elseif l == (n + 1) * 2^nr
        return 2
    else
        throw("wrong length")
    end
end

"""
KL(x,y)

Kullback-Leibler distance
"""
KL(x, y) = sum(x .* (log.(max.(x, eps(Float64))) - log.(max.(y, eps(Float64)))))

"""
trim_hist(h::Array,nh::Array)

"""
function trim_hist(h::Array, nh::Array)
    for i in eachindex(h)
        h[i] = h[i][1:nh[i]]
    end
    return h
end

"""
    distribution_array(param::Vector, cv, dist=Normal)

"""
function distribution_array(param::Vector, cv, dist=Normal)
    d = []
    for i in eachindex(param)
        if dist == LogNormal
            push!(d, LogNormal_meancv(param[i], cv[i]))
        elseif dist == Gamma
            push!(d, Gamma_meancv(param[i], cv[i]))
        else
            push!(d, dist(param[i], cv[i]))
        end

    end
    return d
end


"""
LogNormal_array(param,cv)
Prior distribution arrays
"""
LogNormal_array(param, cv) = distribution_array(param, cv, LogNormal)

"""
Gamma_array(param,cv)



"""
Gamma_array(param, cv) = distribution_array(param, cv, Gamma)

"""
function LogNormalBeta_array(param,cv,ind)


"""
function LogNormalBeta_array(param, cv, ind)
    if ind == 1
        d = [setBeta(param[ind], cv[ind])]
    else
        barind = 1:ind-1
        d = LogNormalarray(param[barind], cv[barind])
        push!(d, setBeta(param[ind], cv[ind]))
    end
    return d
end

"""
    distributionBeta_array(param::Vector, cv::Vector, ind::Int, dist=LogNormal)

fill an array with dist(param,cv)
"""
function distributionBeta_array(param::Vector, cv::Vector, ind::Int, dist=LogNormal)
    if ind == 1
        d = [Beta_meancv(param[ind], cv[ind])]
    else
        barind = 1:ind-1
        d = distribution_array(param[barind], cv[barind], dist)
        push!(d, Beta_meancv(param[ind], cv[ind]))
    end
    return d
end



"""
Gamma_meancv(mean,cv)
LogNormal_meancv(mean,cv)
Beta_meancv(mean,cv)

Reparameterized distributions to mean and coefficient of variation arguments
"""
Gamma_meancv(mean, cv) = Gamma(1 / cv^2, mean * cv^2)
LogNormal_meancv(mean, cv) = LogNormal(mulognormal(mean, cv), sigmalognormal(cv))
function Beta_meancv(m, cv)
    cv2 = cv^2
    fac = (1 - m) / cv2 / m - 1
    if fac <= 0.0
        fac = 2 / m
    end
    alpha = m * fac
    beta = (1 - m) * fac
    Beta(alpha, beta)
end
"""
sigmalognormal(cv)
mulognormal(mean,cv)

compute LogNormal mu and sigma parameters, 
i.e. mean and sigma of exp(N(mu,sigma))
given desired mean and coefficient of variation (cv)
"""
sigmalognormal(cv) = sqrt.(log.(1 .+ cv .^ 2))

mulognormal(mean, cv) = log.(mean) - 0.5 * log.(1 .+ cv .^ 2)

mean_dwelltime(x, t) = t' * x

"""
mean_histogram(x)

"""
mean_histogram(x::Vector{Float64}) = collect(0:length(x)-1)' * x / sum(x)

function mean_histogram(x::Vector{Array})
    y = Vector{Float64}(undef, length(x))
    for i in eachindex(x)
        y[i] = mean_histogram(x[i])
    end
    return y
end

"""
    Difference_Zscore(x1,x2,sig1,sig2) = (x1 .- x2)./ sqrt(sig1.^2 + sig2.^2)


"""
Difference_Zscore(x1, x2, sig1, sig2) = (x1 - x2) / sqrt(sig1^2 + sig2^2)

"""
m2_histogram(x)

"""
m2_histogram(x) = (collect(0:length(x)-1) .^ 2)' * x / sum(x)

"""
var_histogram(x)

"""
function var_histogram(x)
    m2_histogram(x) - mean_histogram(x)^2
end

function moment_histogram(x, n)
    y = collect(0:length(x)-1)
    v = (y .- mean_histogram(x)) .^ n
    v' * x / sum(x)
end

function factorial_moment(h::Vector, n)
    m = 0
    for i in n:length(h)
        a = i - 1
        for j in 1:n-1
            a *= i - j - 1
        end
        m += h[i] * a
    end
    return m / sum(h)
end

function moment_param_estimates(h)
    e1 = factorial_moment(h, 1)
    e2 = factorial_moment(h, 2)
    e3 = factorial_moment(h, 3)
    r1 = e1
    r2 = e2 / e1
    r3 = e3 / e2
    println(r1, ", ", r2, ", ", r3)

    lambda = 2 * r1 * (r3 - r2) / (r1 * r2 - 2 * r1 * r3 + r2 * r3)
    mu = 2 * (r2 - r1) * (r1 - r3) * (r3 - r2) / (r1 * r2 - 2 * r1 * r3 + r2 * r3) / (r1 - 2 * r2 + r3)
    nu = (-r1 * r2 + 2 * r1 * r3 - r2 * r3) / (r1 - 2 * r2 + r3)

    return lambda, mu, nu
end
"""
tstat_2sample(x1,x2)

Compute t statistics of histograms x1 and x2
"""
tstat_2sample(x1, x2) = (mean_histogram(x1) - mean_histogram(x2)) / (sqrt(var_histogram(x1) / length(x1) + var_histogram(x2) / length(x2)))


"""
log_2sample(x1,x2)

compute the log of the ratio of means of histograms x1 and x2
"""
log_2sample(x1, x2) = log(mean_histogram(x1)) - log(mean_histogram(x2))

"""
delta_2sample(x1,x2)

compute the difference of means of x1 and x2
"""
delta_2sample(x1, x2) = mean_histogram(x1) - mean_histogram(x2)

"""
delta_2frac(x1,x2)

compute E(x1) - E(x2) of means normalize dot mean of x1
"""
delta_2frac(x1, x2) = delta_2sample(x1, x2) / mean_histogram(x1)


"""
mediansmooth(xin,window)

median smooth histogram xin over a window
"""
function mediansmooth(xin, window)
    x = copy(xin)
    halfwin = div(window - 1, 2)
    for i in 1:length(x)-window+1
        x[i+halfwin] = median(x[i:i+window-1])
    end
    x
end


"""
residenceprob_G(file,G,header)
residenceprob_G(r,n)

Residence probability of G states
given by exact steady state solution
of master equation

"""
function residenceprob_G(file::String, G, header=false)
    r = get_all_rates(file, header)
    m = size(r)[1]
    p = Array{Any,2}(undef, m, G + 1)
    n = G - 1
    for i in 1:m
        p[i, 1] = r[i, 1]
        p[i, 2:end] = residenceprob_G(r[i, 2:2*n+1], n)
    end
    return p
end

function residenceprob_G(r::Vector, n::Int)
    Gss = Array{Float64,2}(undef, 1, n + 1)
    Gss[1, 1] = 1.0
    for k in 1:n
        Gss[1, k+1] = Gss[1, k] * r[2*k-1] / r[2*k]
    end
    Gss ./= sum(Gss)
end

"""
splicesiteusage()

splice site usage probability is the probabilty of ejection times
the probability of not ejection prior to that point
"""
splicesiteusage(model::GRSMmodel) = splicesiteusage(model.rates, model.G - 1, model.R)
function splicesiteusage(r::Vector, n::Int, nr::Int)
    nu = get_nu(r, n, nr)
    eta = get_eta(r, n, nr)
    ssf = zeros(nr)
    survival = 1
    for j in 1:nr
        ssf[j] = eta[j] / (nu[j+1] + eta[j]) * survival
        survival *= nu[j+1] / (nu[j+1] + eta[j])
    end
    return ssf
end

"""
    onstate_prob(r,components::TComponents)

"""
function onstate_prob(r, components::TComponents, reporter_per_state)
    Qtr = make_mat(components.elementsT, r, components.nT)
    p0 = normalized_nullspace(Qtr)
    return sum(p0[reporter_per_state.>0]), p0
end

onstate_prob(param, model::AbstractGmodel) = onstate_prob(get_rates(param, model), model.components, model.reporter.per_state)

onstate_prob(model::AbstractGmodel) = onstate_prob(model.rates, model.components, model.reporter.per_state)


"""
burstoccupancy(n,nr,r)
Burst size distribution of GRS  model
for total pre-mRNA occupancy and
unspliced (visible) occupancy
obtained by marginalizing over conditional steady state distribution
"""
burstoccupancy(model::GRSMmodel) = burstoccupancy(model.G - 1, model.R, model.rates)

function burstoccupancy(n::Int, nr::Int, r::Vector)
    T = mat_GSR_T(r, n, nr)
    pss = normalized_nullspace(T)
    Rss = zeros(nr)
    Rssvisible = zeros(nr)
    ssf = zeros(nr)
    asum = 0
    for w = 1:nr
        for i = 1:n+1, z = 1:3^nr
            zdigits = digit_vector(z, 3, nr)
            a = i + (n + 1) * (z - 1)
            if sum(zdigits .== 2) == w
                Rssvisible[w] += pss[a]
            end
            if sum(zdigits .> 0) == w
                Rss[w] += pss[a]
            end
        end
    end
    Rss ./= sum(Rss), Rssvisible ./= sum(Rssvisible)
end

model2_mean(ron, roff, eject, decay, nalleles) = 2 * model2_mean(ron, roff, eject, decay)

model2_mean(ron, roff, eject, decay) = ron / (ron + roff) * eject / decay

model2_variance(ron, roff, eject, decay, nalleles) = 2 * model2_variance(ron, roff, eject, decay)

model2_variance(ron, roff, eject, decay) = ron / (ron + roff) * eject / decay + ron * roff / (ron + roff)^2 * eject^2 / decay / (ron + roff + decay)

"""
test_steadystatemodel(model,nhist)

Compare chemical master solution to Gillespie simulation for steadystate mRNA distribution

"""
function test_steadystatemodel(model::AbstractGMmodel, nhist)
    G = model.G
    r = model.rates
    g1 = steady_state(r[1:2*G], G - 1, nhist, model.nalleles)
    g2 = simulatorGM(r[1:2*G], G - 1, nhist, model.nalleles)
    return g1, g2
end

function test_model(data::RNAOnOffData, model::GRSMmodel)
    telegraphsplice0(data.bins, data.nRNA, model.G - 1, model.R, model.rates, 1000000000, 1e-5, model.nalleles)
end
"""
 make_histograms(folder,file,label)
 make_histogram(r)


"""
function make_histograms(folder, file, label)
    if ~ispath(folder)
        mkpath(folder)
    end
    a, h = readdlm(file, ',', header=true)
    for r in eachrow(a)
        f = open("$folder/$(r[1])_$label.txt", "w")
        a = r[2:end]
        a = a[(a.!="").&(a.!="NA")]
        h = make_histogram(Int.(a) .+ 1)
        writedlm(f, h)
        close(f)
    end
end

function make_histogram(r)
    nhist = maximum(r)
    h = zeros(Int, nhist + 1)
    for c in r
        h[c] += 1
    end
    h
end


"""
    expv(v::Array)

exponential of array
"""
expv(v::Array) = exp.(v)


"""
    logv(v::Array)

TBW
"""
logv(v::Array) = log.(v)


"""
    log_shift(v::Float64,a::Float64)
    log_shift(v::Array,a::Float64)

TBW
"""
log_shift(v::Float64, a::Float64) = log(v + a)

log_shift(v::Array, a::Float64) = log.(v .+ a)

"""
    log_shift1(v)

TBW
"""
log_shift1(v) = log_shift(v, 1.0)


"""
    invlog_shift1(v::Float64)

TBW
"""
invlog_shift1(v::Float64) = exp(v) - 1
invlog_shift1(v::Array) = exp.(v) .- 1

"""
    logit(x::Float64)

logit transform
"""
logit(x::Float64) = log(x) - log(1 - x)
logit(x::Array) = log.(x) .- log.(1 .- x)

"""
    invlogit(x::Float64)

inverse logit transform
"""
invlogit(x::Float64) = 1 / (1 + exp(-x))
invlogit(x::Array) = 1 ./ (1 .+ exp.(-x))
"""
logsumexp(u,v)

returns log of the sum of exponentials of u and v

"""
function logsumexp(u::Array, v::Array)
    w = max(u, v)
    if w == -Inf
        return -Inf
    else
        return w .+ log.(exp.(u - w) + exp.(v - w))
    end
end

function logsumexp(u::Float64, v::Float64)
    w = max(u, v)
    if w == -Inf
        return -Inf
    else
        return w + log(exp(u - w) + exp(v - w))
    end
end

"""
logsumexp(v)

returns log of the sum of exponentials of elements of v

"""
function logsumexp(v::Array)
    w = maximum(v)
    if w == -Inf
        return -Inf
    else
        return w .+ log(sum(exp.(v .- w)))
    end
end

function mean_elongationtime(r, transitions, R)
    if R > 0
        n = length(transitions)
        return meantime(r[n+2:n+R+1])
    else
        return 1.0
    end
end

function meantime(r::Vector)
    sum(1 ./ r)
end

"""
    make_array(v)

convert containuer of arrays into a single array
"""
function make_array(v)
    vconcat = Float64[]
    for v in v
        vconcat = [vconcat; v]
    end
    return vconcat
end

"""
    makestring(v::Vector)

convert a vector of strings into a single string
"""
function makestring(v)
    s = ""
    for i in v
        s *= i
    end
    return s
end


#
# function fit_rna_test(root)
#     gene = "CENPL"
#     cell = "HCT116"
#     fish = false
#     data = data_rna(gene,"MOCK","data/HCT116_testdata",fish,"scRNA_test",root)
#     model = model_rna(gene,cell,2,fish,.01,[1,2,3],(),"scRNA_test","scRNA_test",1,".",data,.05,1.0)
#     options = MHOptions(10000,2000,0,120.,1.,100.)
#     fit,stats,waic = run_mh(data,model,options,1);
#     return stats.meanparam, fit.llml
# end
#

# function online_covariance(data1, data2)
#     meanx = meany = C = n = 0
#     for x in data1, y in data2
#         n += 1
#         dx = x - meanx
#         meanx += dx / n
#         meany += (y - meany) / n
#         C += dx * (y - meany)
#         population_covar = C / n
#         # Bessel's correction for sample variance
#         sample_covar = C / (n - 1)
#     end
# end
