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

pooled_mean(means,counts) = counts'*means/sum(counts)
pooled_variance(vars,counts) = (counts .- 1)'*vars/sum(counts .- 1)
pooled_std(std,counts) = sqrt(pooled_variance(std.^2,counts))

"""
ternary(x)
convert digits of base 3 number to base 10
"""
function ternary(x::Vector)
    nr = length(x)
    x' * 3 .^ collect(0:nr-1)
end
"""
sigmalognormal(cv)
compute LogNormal sigma parameter given desired coefficient of variation
"""
sigmalognormal(cv) = sqrt.(log.(1 .+ cv .^ 2))

KL(x,y) = sum(x .* (log.(max.(x,eps(Float64))) - log.(max.(y,eps(Float64)))))

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
