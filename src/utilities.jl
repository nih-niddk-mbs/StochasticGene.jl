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

function findbase(l,n,nr)
    if l == (n+1)*3^nr
        return 3
    elseif l == (n+1)*2^nr
        return 2
    else
        throw("wrong length")
    end
end

KL(x,y) = sum(x .* (log.(max.(x,eps(Float64))) - log.(max.(y,eps(Float64)))))


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
Parameter information
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
