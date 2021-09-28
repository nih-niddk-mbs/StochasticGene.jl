"""
Strings
"""
const rchar = "rml"
const rmchar = "rmean"
const rqchar = "rquant"
const wchar = "waic"
const routchar = "rout"
const runchar = "run"
const chichar = "chi"
const chioutchar = "chiout"
const rlastchar = "rlast"
const onchar = "LCON"
const offchar = "LCOFF"
const fishchar = "FISH"
const txtstr = ".txt"

"""
eig_decompose(M)
"""
function eig_decompose(M)
    Meiv = eigen(M)
    return Meiv.values, Meiv.vectors
end

"""
nonzero_rows(T)
Locate all rows that have at least one nonzero element
"""
function nonzero_rows(T)
    n = size(T)[2]
    nonzeros = Array{Int}(undef,0)
    for  a in 1:size(T)[1]
        if T[a,:] != zeros(n)
            push!(nonzeros,a)
        end
    end
    return nonzeros
end

normalize_histogram(hist) = hist/sum(hist)


"""
var_update(vartuple, newValue)
iteration step for online algorithm to compute
mean and variance
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
"""
function cov_update(covtuple,newValue1,newValue2)
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
truncatehistogram(x::Array{Int,1},yield::Float64,nhistmax::Int)
"""
function truncate_histogram(x::Array,yield::Float64,nhistmax::Int)
    counts = sum(x)
    if counts > 0
        nhist = 1
        ratio = 0.
        # Only keep yield of all data to limit support of histogram
        while ratio < yield && nhist <= nhistmax
            ratio = sum(x[1:nhist])/counts
            nhist += 1
        end
        return x[1:min(nhist,length(x))]
    else
        return 0
    end
end

function combine_histogram(x1::Array,x2::Array)
    n = min(length(x1),length(x2))
    x1[1:n] + x2[1:n]
end

pooled_mean(means,counts) = counts'*means/sum(counts)
pooled_variance(vars,counts) = (counts .- 1)'*vars/sum(counts .- 1)
pooled_std(std,counts) = sqrt(pooled_variance(std.^2,counts))

"""
decimal(x::Vector)
convert digits of base 3 number to base 10
"""
function decimal(x::Vector,base::Int=3)
    nr = length(x)
    x' * base .^ collect(0:nr-1)
end
"""
findbase(l,n,nr)
Find the number of G states for transition rate matrix of size l
"""
function findbase(l,n,nr)
    if l == (n+1)*3^nr
        return 3
    elseif l == (n+1)*2^nr
        return 2
    else
        throw("wrong length")
    end
end

"""
KL(x,y)

Kullback-Leibler distance
"""
KL(x,y) = sum(x .* (log.(max.(x,eps(Float64))) - log.(max.(y,eps(Float64)))))

"""
trim(h::Array,nh::Array)

"""
function trim(h::Array,nh::Array)
    for i in eachindex(h)
        h[i] = h[i][1:nh[i]]
    end
    return h
end

"""
LogNormal_array(param,cv)
Prior distribution arrays
"""
LogNormal_array(param,cv) = distributionarray(param,cv,LogNormal)
Gamma_array(param,cv) = distributionarray(param,cv,Gamma)

function LogNormalBeta_array(param,cv,ind)
    if ind == 1
        d = [setBeta(param[ind],cv[ind])]
    else
        barind = 1:ind-1
        d = LogNormalarray(param[barind],cv[barind])
        push!(d,setBeta(param[ind],cv[ind]))
    end
    return d
end
"""
distribution_array(param,cv,dist)

fill an array with dist(param,cv)
"""

function distributionBeta_array(param::Vector,cv::Vector,ind::Int,dist=LogNormal)
    if ind == 1
        d = [Beta_meancv(param[ind],cv[ind])]
    else
        barind = 1:ind-1
        d = distribution_array(param[barind],cv[barind],dist)
        push!(d,Beta_meancv(param[ind],cv[ind]))
    end
    return d
end
function distribution_array(param::Vector,cv,dist=LogNormal)
    d = []
    for i in eachindex(param)
        push!(d,dist(param[i],cv[i]))
    end
    return d
end

"""
Gamma_meancv(mean,cv)
LogNormal_meancv(mean,cv)
Beta_meancv(mean,cv)

Reparameterized distributions to mean and coefficient of variation arguments
"""
Gamma_meancv(mean,cv) = Gamma(1/cv^2,mean * cv^2)
LogNormal_meancv(mean,cv) = LogNormal(mulognormal(mean,cv), sigmalognormal(cv))
function Beta_meancv(m,cv)
    cv2 = cv^2
    fac = (1-m)/cv2/m - 1
    if fac <= 0.
        fac = 2/m
    end
    alpha = m*fac
    beta = (1-m)*fac
    Beta(alpha,beta)
end
"""
sigmalognormal(cv)
mulognormal(mean,cv)

compute LogNormal mu and sigma parameters
given desired mean and coefficient of variation
"""
sigmalognormal(cv) = sqrt.(log.(1 .+ cv .^ 2))

mulognormal(mean,cv) = log.(mean) - .5*log.(1 .+ cv .^ 2)

"""
Parameter information for GR models
"""
function num_rates(G,R)
    n = G - 1
    2*n + 2*R + 2
end
function Grange(G)
    n = G - 1
    1 : 2*n
end
function initiation(G)
    n = G - 1
    2*n + 1
end
function Rrange(G,R)
    n = G - 1
    2*n + 2 : 2*n + R + 1
end
function Srange(G,R)
    n = G - 1
    2*n + R + 2 : 2*n + 2*R + 1
end

mean_histogram(x) = (collect(1:length(x)) .- 1)' * x/sum(x)
m2_histogram(x) = ((collect(1:length(x)) .- 1).^2)' * x/sum(x)
function var_histogram(x)
    m2_histogram(x) - mean_histogram(x)^2
end
# mean_histogram(x) = ((collect(length(x)) .- 1)' * x)/sum(x)
tstat_2sample(x1,x2) = (mean_histogram(x1) - mean_histogram(x2))/(sqrt(var_histogram(x1)/length(x1) + var_histogram(x2)/length(x2)))

log_2sample(x1,x2) = log(mean_histogram(x1)) - log(mean_histogram(x2))
delta_2sample(x1,x2) = mean_histogram(x1)-mean_histogram(x2)
delta_2frac(x1,x2) = delta_2sample(x1,x2)/mean_histogram(x1)

function mediansmooth(xin,window)
    x = copy(xin)
    halfwin = div(window-1,2)
    for i in 1:length(x)-window+1
        x[i+halfwin] = median(x[i:i+window-1])
    end
    x
end

"""
plot_histogram()

functions to plot data and model predicted histograms

"""
function plot_histogram(data::RNAData,model::GMlossmodel)
    h=likelihoodarray(model.rates,data,model)
    for i in eachindex(h)
        figure()
        plot(h[i])
        plot(normalize_histogram(data.histRNA[i]))
    end
    return h
end
function plot_histogram(data::AbstractRNAData{Array{Array,1}},model)
    h=likelihoodarray(model.rates,data,model)
    figure(data.gene)
    for i in eachindex(h)
        plot(h[i])
        plot(normalize_histogram(data.histRNA[i]))
    end
    return h
end
function plot_histogram(data::AbstractRNAData{Array{Float64,1}},model)
    h=likelihoodfn(model.rates,data,model)
    figure(data.gene)
    plot(h)
    plot(normalize_histogram(data.histRNA))
    return h
end

function plot_histogram(data::RNALiveCellData,model)
    h=likelihoodtuple(model.rates,data,model)
    figure(data.gene)
    plot(h[1])
    plot(normalize_histogram(data.OFF))
    plot(h[2])
    plot(normalize_histogram(data.ON))
    figure("FISH")
    plot(h[3])
    plot(normalize_histogram(data.histRNA))
    return h
end

function plot_histogram(data::TransientRNAData,model::AbstractGMmodel)
    h=StochasticGene.likelihoodarray(model.rates,data,model)
    for i in eachindex(h)
        figure(data.gene *":T" * "$(data.time[i])")
        plot(h[i])
        plot(normalize_histogram(data.histRNA[i]))
    end
    return h
end

function plot_histogram(data::RNAData,model::AbstractGMmodel)
    h=StochasticGene.likelihoodfn(model.rates,data,model)
    figure(data.gene)
    plot(h)
    plot(normalize_histogram(data.histRNA))
    return h
end

"""
residenceprob_G(file,G,header)
residenceprob_G(r,n)

Residence probability of G states
given by exact steady state solution
of master equation

"""
function residenceprob_G(file::String,G,header=false)
    r = get_all_rates(file,header)
    m = size(r)[1]
    p = Array{Any,2}(undef,m,G+1)
    n = G-1
    for i in 1:m
        p[i,1] = r[i,1]
        p[i,2:end] = residenceprob_G(r[i,2:2*n+1],n)
    end
    return p
end

function residenceprob_G(r::Vector,n::Int)
    Gss = Array{Float64,2}(undef,1,n+1)
    Gss[1,1] = 1.
    for k in 1:n
        Gss[1,k+1] = Gss[1,k]*r[2*k-1]/r[2*k]
    end
    Gss ./= sum(Gss)
end

"""
splicesiteusage()

splice site usage probability is the probabilty of ejection times
the probability of not ejection prior to that point
"""
splicesiteusage(model::GRSMmodel) = splicesiteusage(model.rates,model.G-1,model.R)
function splicesiteusage(r::Vector,n::Int,nr::Int)
    nu = get_nu(r,n,nr)
    eta = get_eta(r,n,nr)
	ssf = zeros(nr)
	survival = 1
	for j in 1:nr
		ssf[j] = eta[j]/(nu[j+1]+eta[j])*survival
		survival *= nu[j+1]/(nu[j+1]+eta[j])
	end
	return ssf
end

"""
burstsize(n,nr,r)
Burst size distribution of GRS  model
for total pre-mRNA occupancy and
unspliced (visible) occupancy
obtained by marginalizing over conditional steady state distribution
"""
burstsize(model::GRSMmodel) = burstsize(model.G-1,model.R,model.rates)
function burstsize(n::Int,nr::Int,r::Vector)
    T =  mat_GSR_T(r,n,nr)
    pss = normalized_nullspace(T)
	Rss = zeros(nr)
	Rssvisible = zeros(nr)
	ssf = zeros(nr)
	asum = 0
	for w=1:nr
		for i = 1:n+1, z = 1:3^nr
			zdigits = digits(z-1,base=3,pad=nr)
			a = i + (n+1)*(z-1)
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

"""
teststeadystatemodel(model,nhist)

Compare chemical master solution to Gillespie simulation for steadystate mRNA distribution

"""
function teststeadystatemodel(model::AbstractGMmodel,nhist)
    G = model.G
    r = model.rates
    g1 = steady_state(r[1:2*G],G-1,nhist,model.nalleles)
    g2 = telegraph(G-1,r[1:2*G],10000000,1e-5,nhist,model.nalleles)
    return g1,g2
end

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
