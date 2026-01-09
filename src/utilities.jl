# This file is part of StochasticGene.jl    

# utilities.jl


"""
    invert_dict(D)

Inverts a dictionary.

# Arguments
- `D`: The dictionary to be inverted.

# Description
This function inverts a dictionary by swapping its keys and values.

# Returns
- `Dict`: The inverted dictionary.
"""
invert_dict(D) = Dict(D[k] => k for k in keys(D))



"""
    eig_decompose(M)

Performs eigenvalue decomposition on a matrix.

# Arguments
- `M`: The matrix to be decomposed.

# Description
This function takes a matrix `M` and returns its eigenvalues and eigenvectors.

# Returns
- `Tuple{Vector{Float64}, Matrix{Float64}}`: A tuple containing the eigenvalues and eigenvectors of the matrix.
"""
function eig_decompose(M)
    Meiv = eigen(M)
    return Meiv.values, Meiv.vectors
end

"""
    crosscorrelation_function(x, y, lags; meanx=0.0, meany=0.0)

Explicitly compute the unbiased cross-correlation function without calling Julia's StatsBase.crosscov.

The cross-correlation function is defined as:
    R_XY(τ) = (1/(T-τ)) * Σ_{t=1}^{T-τ} (X(t) - μ_X)(Y(t+τ) - μ_Y)

If `meanx=0.0` and `meany=0.0` (default), computes the uncentered cross-correlation:
    R_XY(τ) = (1/(T-τ)) * Σ_{t=1}^{T-τ} X(t)Y(t+τ)

The normalization `1/(T-τ)` is used for unbiased estimation, matching the number of valid pairs at each lag.

# Arguments
- `x`: First time series vector
- `y`: Second time series vector  
- `lags`: Vector of time lags (can include negative values)
- `meanx::Float64=0.0`: Mean to subtract from x (if 0.0, no centering)
- `meany::Float64=0.0`: Mean to subtract from y (if 0.0, no centering)

# Returns
- `Vector{Float64}`: Unbiased cross-correlation function at specified lags

# Example
```julia
x = randn(100)
y = randn(100)
lags = collect(-20:20)
# Uncentered cross-correlation
R_XY = crosscorrelation_function(x, y, lags)
# Cross-covariance with empirical means
R_XY_centered = crosscorrelation_function(x, y, lags, meanx=mean(x), meany=mean(y))
# Cross-covariance with theoretical means
R_XY_theory = crosscorrelation_function(x, y, lags, meanx=μ_X_theory, meany=μ_Y_theory)
```
"""
function crosscorrelation_function(x, y, lags; meanx::Float64=0.0, meany::Float64=0.0)
    n = length(x)
    if length(y) != n
        error("x and y must have the same length. Got length(x)=$n, length(y)=$(length(y))")
    end
    
    # Center the data if means are provided
    if meanx != 0.0 || meany != 0.0
        x = x .- meanx
        y = y .- meany
    end
    
    # Pre-allocate result
    result = Vector{Float64}(undef, length(lags))
    
    # Compute cross-correlation for each lag
    for (i, τ) in enumerate(lags)
        τ_abs = abs(τ)
        n_valid = n - τ_abs
        
        if n_valid <= 0
            result[i] = 0.0
            continue
        end
        
        # Compute sum of products for valid pairs
        sum_xy = 0.0
        if τ >= 0
            # Positive lag: X(t) * Y(t+τ)
            for t in 1:n_valid
                sum_xy += x[t] * y[t + τ]
            end
        else
            # Negative lag: X(t-|τ|) * Y(t) = Y(t) * X(t-|τ|)
            for t in 1:n_valid
                sum_xy += y[t] * x[t + τ_abs]
            end
        end
        
        # Normalize by number of valid pairs: 1/(T-τ)
        result[i] = sum_xy / n_valid
    end
    
    return result
end

"""
    unbiased_crosscov(x, y, lags; demean=true)

Compute unbiased uncentered cross-correlation or cross-covariance that accounts for different numbers of valid pairs at different lags.

This function wraps `StatsBase.crosscov` and applies a correction factor to account for the fact
that `StatsBase.crosscov` normalizes by `n` (length of input) for all lags, but the actual number
of valid pairs for lag τ is `n - |τ|`. The correction factor is `n / (n - |τ|)`.

When `demean=false`, computes the uncentered cross-correlation function:
    R_XY(τ) = (1/(T-τ)) * Σ_{t=1}^{T-τ} X(t)Y(t+τ)

When `demean=true`, computes the cross-covariance function:
    C_XY(τ) = (1/(T-τ)) * Σ_{t=1}^{T-τ} (X(t) - μ_X)(Y(t+τ) - μ_Y)

where μ_X and μ_Y are computed over the full time series (not windowed means).

The normalization `1/(T-τ)` is used for unbiased estimation, matching the number of valid pairs at each lag.

# Arguments
- `x`: First time series vector
- `y`: Second time series vector  
- `lags`: Vector of time lags (can include negative values)
- `demean`: Whether to remove means before computation (default: true). If false, returns uncentered cross-correlation R_XY(τ).

# Returns
- `Vector{Float64}`: Unbiased uncentered cross-correlation (if `demean=false`) or cross-covariance (if `demean=true`) at specified lags

# Example
```julia
x = randn(100)
y = randn(100)
lags = collect(-20:20)
# Uncentered cross-correlation
R_XY = unbiased_crosscov(x, y, lags, demean=false)
# Cross-covariance
C_XY = unbiased_crosscov(x, y, lags, demean=true)
```
"""
function unbiased_crosscov(x, y, lags; demean::Bool=true)
    # Compute empirical means if demeaning is requested
    if demean
        meanx = mean(x)
        meany = mean(y)
    else
        meanx = 0.0
        meany = 0.0
    end
    
    # Call crosscorrelation_function with computed means
    return crosscorrelation_function(x, y, lags; meanx=meanx, meany=meany)
end

"""
    nonzero_rows(T)

Returns an array of row indices that have at least one nonzero element.

# Arguments
- `T`: The matrix to be checked.

# Description
This function returns an array of row indices that have at least one nonzero element for the given matrix `T`.

# Returns
- `Vector{Int}`: An array of row indices with at least one nonzero element.
"""
function nonzero_rows(T::AbstractMatrix)
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
    nonzero_states(states, nonzeros)

Returns a reindexed state vector with zeros removed.

# Arguments
- `states`: The original state vector.
- `nonzeros`: The vector of nonzero indices.

# Description
This function returns a reindexed state vector with zeros removed by finding the intersection of `states` and `nonzeros`.

# Returns
- `Vector{Int}`: The reindexed state vector with zeros removed.
"""
nonzero_states(states, nonzeros) = [findfirst(o .== nonzeros) for o in intersect(states, nonzeros)]

"""
    normalize_histogram(hist)

Returns a normalized histogram.

# Arguments
- `hist`: The histogram to be normalized.

# Description
This function returns a normalized version of the given histogram `hist`.

# Returns
- `Vector{Float64}`: The normalized histogram.
"""
normalize_histogram(hist) = hist / sum(hist)


"""
    var_update(vartuple::Tuple, newValue)

Updates the variance calculation with a new value using Welford's online algorithm.

# Arguments
- `vartuple`: A tuple containing the current count, mean, and M2 (sum of squares of differences from the current mean).
- `newValue`: The new value to be added to the variance calculation.

# Description
This function updates the variance calculation with a new value using Welford's online algorithm. It returns the updated count, mean, and M2.

# Returns
- `Tuple{Int, Float64, Float64}`: The updated count, mean, and M2.
"""
function var_update(vartuple::Tuple, newValue)
    count, mean, M2 = vartuple
    count += 1
    delta = newValue - mean
    mean += delta / count
    delta2 = newValue - mean
    M2 += delta * delta2
    return count, mean, M2
end

function var_update(vartuple::Tuple, newValue::Vector{Float64})
    count, mean, M2 = vartuple
    count += 1
    delta = newValue .- mean
    mean += delta ./ count
    delta2 = newValue .- mean
    M2 += delta .* delta2
    return count, mean, M2
end

function test_varupdate(n)
    a = (0, zeros(3), zeros(3))
    for i in 1:n
        a = var_update(a, randn(3))
    end
    count, mean, M2 = a
    return mean, M2 / (count - 1)
end

"""
    cov_update(covtuple, newValue1, newValue2)

Updates the covariance calculation with new values using an online algorithm.

# Arguments
- `covtuple`: A tuple containing the current count, mean of x, mean of y, and C (sum of products of differences from the current means).
- `newValue1`: The new value for the first variable.
- `newValue2`: The new value for the second variable.

# Description
This function updates the covariance calculation with new values using an online algorithm. It returns the updated count, mean of x, mean of y, and C.

# Returns
- `Tuple{Int, Float64, Float64, Float64}`: The updated count, mean of x, mean of y, and C.
"""
function cov_update(covtuple::Tuple, newValue1, newValue2)
    count, meanx, meany, C = covtuple
    count += 1
    dx = newValue1 - meanx
    meanx .+= dx / count
    meany .+= (newValue2 - meany) / count
    C .+= dx * (newValue2 - meany)
    return count, meanx, meany, C
end

"""
    mean_update(meantuple::Tuple, newValue)

Updates the mean with a new value using Welford's online algorithm.

# Arguments
- `meantuple`: A tuple containing the current mean and count.
- `newValue`: The new value to be added to the mean calculation.

# Description
This function updates the mean with a new value using Welford's online algorithm. It returns the updated mean and count.

# Returns
- `Tuple{Float64, Int}`: The updated mean and count.
"""
function mean_update(meantuple::Tuple, newValue)
    mean, count = meantuple
    count += 1
    delta = newValue - mean
    mean += delta / count
    return mean, count
end


"""
    mean_update(mean::Vector{Float64}, count::Int, newValue::Vector{Float64})

Updates the mean vector with a new value using Welford's online algorithm.

# Arguments
- `mean`: The current mean vector.
- `count`: The current count of values.
- `newValue`: The new value vector to be added to the mean calculation.

# Description
This function updates the mean vector with a new value using Welford's online algorithm. It returns the updated mean vector and count.

# Returns
- `Tuple{Vector{Float64}, Int}`: The updated mean vector and count.
"""
function mean_update(mean::Vector{Float64}, count::Int, newValue::Vector{Float64})
    delta = newValue .- mean
    mean .+= delta / count
    return mean
end

function mean_update(mean::Float64, count, newValue)
    delta = newValue - mean
    mean += delta / count
    return mean
end


"""
    finalize(vartuple)

Retrieves the mean, variance, and sample variance from an aggregate.

# Arguments
- `vartuple`: A tuple containing the count, mean, and M2 (sum of squares of differences from the current mean).

# Description
This function retrieves the mean, variance, and sample variance from an aggregate. It uses the count, mean, and M2 to compute the variance and sample variance.

# Returns
- `Tuple{Float64, Float64, Float64}`: The mean, variance, and sample variance.
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

Returns a truncated histogram that accounts for a specified fraction of the total sum.

# Arguments
- `x`: The histogram array.
- `yield`: The fraction of the total sum to retain.
- `nhistmax`: The maximum number of histogram bins to retain.

# Description
This function returns a truncated version of the histogram `x` that accounts for the specified fraction `yield` of the total sum. It iterates through the histogram bins and retains bins until the specified fraction of the total sum is reached or the maximum number of bins `nhistmax` is reached.

# Returns
- `Array`: The truncated histogram.
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
    trim_hist(hist::Array, nRNA::Array)

Trims the histogram to match the length of the RNA data.

# Arguments
- `hist`: The histogram array.
- `nRNA`: The RNA data array.

# Description
This function trims the histogram `hist` to match the length of the RNA data `nRNA`. It ensures that the histogram and RNA data arrays have the same length.

# Returns
- `Array`: The trimmed histogram.
"""
function trim_hist(h::Array, nh::Array)
    for i in eachindex(h)
        h[i] = h[i][1:nh[i]]
    end
    return h
end


"""
    combine_histogram(x1::Array, x2::Array)

Adds two histograms element-wise.

# Arguments
- `x1`: The first histogram array.
- `x2`: The second histogram array.

# Description
This function adds two histograms `x1` and `x2` element-wise. It returns a new histogram that is the element-wise sum of the two input histograms. The length of the resulting histogram is the minimum of the lengths of the input histograms.

# Returns
- `Array`: The combined histogram.
"""
function combine_histogram(x1::Array, x2::Array)
    n = min(length(x1), length(x2))
    x1[1:n] + x2[1:n]
end

"""
    pooled_mean(means, counts)

Computes the weighted average of means given count totals.

# Arguments
- `means`: An array of means.
- `counts`: An array of count totals corresponding to the means.

# Description
This function computes the weighted average of the means given the count totals. The weights are the count totals.

# Returns
- `Float64`: The weighted average of the means.
"""
pooled_mean(means, counts) = counts' * means / sum(counts)

"""
    pooled_variance(vars, counts)

Computes the weighted average of variances given count totals.

# Arguments
- `vars`: An array of variances.
- `counts`: An array of count totals corresponding to the variances.

# Description
This function computes the weighted average of the variances given the count totals. The weights are the count totals.

# Returns
- `Float64`: The weighted average of the variances.
"""
pooled_variance(vars, counts) = (counts .- 1)' * vars / sum(counts .- 1)

"""
    pooled_std(std, counts)

Computes the weighted average of standard deviations given count totals.

# Arguments
- `std`: An array of standard deviations.
- `counts`: An array of count totals corresponding to the standard deviations.

# Description
This function computes the weighted average of the standard deviations given the count totals. The weights are the count totals.

# Returns
- `Float64`: The weighted average of the standard deviations.
"""
pooled_std(std, counts) = sqrt(pooled_variance(std .^ 2, counts))

"""
    var_ratio(mua, mub, vara, varb, cov)

Computes the ratio of variances of two variables given their means, variances, and covariance.

# Arguments
- `mua`: The mean of the first variable.
- `mub`: The mean of the second variable.
- `vara`: The variance of the first variable.
- `varb`: The variance of the second variable.
- `cov`: The covariance of the two variables.

# Description
This function computes the ratio of variances of two variables given their means, variances, and covariance.

# Returns
- `Float64`: The ratio of variances of the two variables.
"""
var_ratio(mua, mub, vara, varb, cov) = mua^2 / mub^2 * (vara / mua^2 - 2 * cov / (mua * mub) + varb / mub^2)

"""
    decimal(x::Vector, base::Int=3)

Converts a vector of digits in a specified base to a decimal (base 10) number.

# Arguments
- `x`: A vector representing the number in the original base, with the least significant digit first.
- `base`: The base of the original number (default is 3).

# Description
This function converts a number represented as a vector `x` from the specified base to a decimal (base 10) number. The vector `x` should have the least significant digit first.

# Returns
- `Int`: The number converted to base 10.
"""
function decimal(x::Vector, base::Int=3)
    nr = length(x)
    x' * base .^ collect(0:nr-1)
end

"""
    digit_vector(z, base, R)

Converts a number to a vector of digits in a specified base.

# Arguments
- `z`: The number to be converted.
- `base`: The base for the conversion.
- `R`: The length of the resulting vector.

# Description
This function converts the number `z` to a vector of length `R` in the specified `base`.

# Returns
- `Vector{Int}`: The vector of digits representing the number in the specified base.
"""
digit_vector(z, base, R) = digits(z - 1, base=base, pad=R)

"""
    findbase(l, n, nr)

Finds the number of G states for a transition rate matrix of a given size.

# Arguments
- `l`: The size of the transition rate matrix.
- `n`: The number of states.
- `nr`: The number of rates.

# Description
This function finds the number of G states for a transition rate matrix of size `l`.

# Returns
- `Int`: The number of G states.
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
    KL(x, y)

Computes the Kullback-Leibler (KL) divergence between two probability distributions.

# Arguments
- `x`: The first probability distribution (vector).
- `y`: The second probability distribution (vector).

# Description
This function computes the Kullback-Leibler (KL) divergence between two probability distributions `x` and `y`. The KL divergence is a measure of how one probability distribution diverges from a second, expected probability distribution. The function ensures that the values in `x` and `y` are not zero by replacing them with a small epsilon value before computing the logarithms.

# Returns
- `Float64`: The KL divergence between the two probability distributions.
"""
KL(x, y) = sum(x .* (log.(max.(x, eps(Float64))) - log.(max.(y, eps(Float64)))))

"""
    distribution_array(param::Vector, cv, dist=Normal)

Creates an array of distributions based on the provided parameters and coefficient of variation (cv).

# Arguments
- `param`: A vector of parameters (means) for the distributions.
- `cv`: A vector of coefficients of variation for the distributions.
- `dist`: The type of distribution to create (default is `Normal`). Supported distributions include `LogNormal`, `Gamma`, and any other distribution that can be constructed with the provided parameters and cv.

# Description
This function creates an array of distributions based on the provided parameters and coefficient of variation (cv). For each parameter in `param`, it creates a distribution of the specified type (`dist`) with the corresponding coefficient of variation from `cv`. If `dist` is `LogNormal` or `Gamma`, it uses specialized functions (`LogNormal_meancv` and `Gamma_meancv`) to create the distributions.

# Returns
- `Vector{Distribution}`: An array of distributions created based on the provided parameters and coefficient of variation.
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

function truncated_normal(μ, σ, k=4)
    lower = μ - k * σ
    upper = μ + k * σ
    truncated(Normal(μ, σ), lower, upper)
end

function sigmanormal(mean, cv)
    return max(abs(mean * cv), abs(cv))
end

"""
    sigmalognormal(cv)

Calculates the standard deviation parameter (σ) for a LogNormal distribution based on the coefficient of variation (cv).

# Arguments
- `cv`: The coefficient of variation of the LogNormal distribution.

# Description
This function calculates the standard deviation parameter (σ) for a LogNormal distribution using the provided coefficient of variation (cv).

# Returns
- `Float64`: The standard deviation parameter (σ) for the LogNormal distribution.
"""
sigmalognormal(cv) = sqrt.(log.(1 .+ cv .^ 2))

"""
    mulognormal(mean, cv)

Calculates the mean parameter (μ) for a LogNormal distribution based on the provided mean and coefficient of variation (cv).

# Arguments
- `mean`: The mean of the LogNormal distribution.
- `cv`: The coefficient of variation of the LogNormal distribution.

# Description
This function calculates the mean parameter (μ) for a LogNormal distribution using the provided mean and coefficient of variation (cv).

# Returns
- `Float64`: The mean parameter (μ) for the LogNormal distribution.
"""
mulognormal(mean, cv) = log.(mean) - 0.5 * log.(1 .+ cv .^ 2)

"""
    Gamma_meancv(mean, cv)

Creates a Gamma distribution based on the provided mean and coefficient of variation (cv).

# Arguments
- `mean`: The mean of the Gamma distribution.
- `cv`: The coefficient of variation of the Gamma distribution.

# Description
This function creates a Gamma distribution based on the provided mean and coefficient of variation (cv). It calculates the shape and scale parameters for the Gamma distribution using the mean and cv.

# Returns
- `Gamma`: A Gamma distribution created based on the provided mean and coefficient of variation.
"""
Gamma_meancv(mean, cv) = Gamma(1 / cv^2, mean * cv^2)

"""
    LogNormal_meancv(mean, cv)

Creates a LogNormal distribution based on the provided mean and coefficient of variation (cv).

# Arguments
- `mean`: The mean of the LogNormal distribution.
- `cv`: The coefficient of variation of the LogNormal distribution.

# Description
This function creates a LogNormal distribution based on the provided mean and coefficient of variation (cv). It calculates the parameters `μ` and `σ` for the LogNormal distribution using the mean and cv.

# Returns
- `LogNormal`: A LogNormal distribution created based on the provided mean and coefficient of variation.
"""
LogNormal_meancv(mean, cv) = LogNormal(mulognormal(mean, cv), sigmalognormal(cv))



"""
    LogNormal_array(param, cv)

Creates an array of LogNormal distributions based on the provided parameters and coefficient of variation (cv).

# Arguments
- `param`: A vector of parameters (means) for the LogNormal distributions.
- `cv`: A vector of coefficients of variation for the LogNormal distributions.

# Description
This function creates an array of LogNormal distributions based on the provided parameters and coefficient of variation (cv). For each parameter in `param`, it creates a LogNormal distribution with the corresponding coefficient of variation from `cv`.

# Returns
- `Vector{LogNormal}`: An array of LogNormal distributions created based on the provided parameters and coefficient of variation.
"""
LogNormal_array(param, cv) = distribution_array(param, cv, LogNormal)

"""
    Gamma_array(param, cv)

Creates an array of Gamma distributions based on the provided parameters and coefficient of variation (cv).

# Arguments
- `param`: A vector of parameters (means) for the Gamma distributions.
- `cv`: A vector of coefficients of variation for the Gamma distributions.

# Description
This function creates an array of Gamma distributions based on the provided parameters and coefficient of variation (cv). For each parameter in `param`, it creates a Gamma distribution with the corresponding coefficient of variation from `cv`.

# Returns
- `Vector{Gamma}`: An array of Gamma distributions created based on the provided parameters and coefficient of variation.
"""
Gamma_array(param, cv) = distribution_array(param, cv, Gamma)

"""
    setBeta(mean, cv)

Creates a Beta distribution based on the provided mean and coefficient of variation (cv).

# Arguments
- `mean`: The mean of the Beta distribution.
- `cv`: The coefficient of variation of the Beta distribution.

# Description
This function creates a Beta distribution based on the provided mean and coefficient of variation (cv). It calculates the parameters `α` and `β` for the Beta distribution using the mean and cv.

# Returns
- `Beta`: A Beta distribution created based on the provided mean and coefficient of variation.
"""
function setBeta(mean, cv)
    var = (cv * mean)^2
    α = mean * (mean * (1 - mean) / var - 1)
    β = (1 - mean) * (mean * (1 - mean) / var - 1)
    return Beta(α, β)
end
"""
    Beta_meancv(mean, cv)

Creates a Beta distribution based on the provided mean and coefficient of variation (cv).

# Arguments
- `mean`: The mean of the Beta distribution.
- `cv`: The coefficient of variation of the Beta distribution.

# Description
This function creates a Beta distribution based on the provided mean and coefficient of variation (cv). It calculates the parameters `α` and `β` for the Beta distribution using the mean and cv.

# Returns
- `Beta`: A Beta distribution created based on the provided mean and coefficient of variation.
"""
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
    LogNormalBeta_array(param, cv, ind)

Creates an array of LogNormal distributions and a Beta distribution based on the provided parameters and coefficient of variation (cv).

# Arguments
- `param`: A vector of parameters (means) for the distributions.
- `cv`: A vector of coefficients of variation for the distributions.
- `ind`: The index at which to insert the Beta distribution.

# Description
This function creates an array of LogNormal distributions and a Beta distribution based on the provided parameters and coefficient of variation (cv). For indices before `ind`, it creates LogNormal distributions. At index `ind`, it creates a Beta distribution using the `setBeta` function.

# Returns
- `Vector{Distribution}`: An array of LogNormal and Beta distributions created based on the provided parameters and coefficient of variation.
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

Creates an array of distributions and a Beta distribution based on the provided parameters and coefficient of variation (cv).

# Arguments
- `param`: A vector of parameters (means) for the distributions.
- `cv`: A vector of coefficients of variation for the distributions.
- `ind`: The index at which to insert the Beta distribution.
- `dist`: The type of distribution to create for indices before `ind` (default is `LogNormal`).

# Description
This function creates an array of distributions and a Beta distribution based on the provided parameters and coefficient of variation (cv). For indices before `ind`, it creates distributions of the specified type (`dist`). At index `ind`, it creates a Beta distribution using the `Beta_meancv` function.

# Returns
- `Vector{Distribution}`: An array of distributions and a Beta distribution created based on the provided parameters and coefficient of variation.
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

# Soft truncated normal: returns a callable logpdf-like function
# Within [μ-kσ, μ+kσ]: logpdf(Normal(μ, σ), x)
# Outside: logpdf at boundary minus a quadratic penalty
function soft_truncated_normal_lpdf(μ, σ, k=4, penalty=10.0)
    lower = μ - k * σ
    upper = μ + k * σ
    norm = Normal(μ, σ)
    function logpdf_soft(x)
        if x < lower
            return logpdf(norm, lower) - penalty * (lower - x)^2
        elseif x > upper
            return logpdf(norm, upper) - penalty * (x - upper)^2
        else
            return logpdf(norm, x)
        end
    end
    return logpdf_soft
end

"""
    mean_dwelltime(x, t)

Calculates the mean dwell time given a vector of states and a vector of times.

# Arguments
- `x`: A vector of states.
- `t`: A vector of times corresponding to the states.

# Description
This function calculates the mean dwell time by taking the dot product of the vector of states `x` and the vector of times `t`.

# Returns
- `Float64`: The mean dwell time.
"""
mean_dwelltime(x, t) = t' * x

"""
    mean_histogram(x::Vector{Float64})

Calculates the mean of a histogram.

# Arguments
- `x`: A vector representing the histogram.

# Description
This function calculates the mean of a histogram `x` by taking the dot product of the vector of bin indices and the histogram values, then dividing by the sum of the histogram values.

# Returns
- `Float64`: The mean of the histogram.
"""
mean_histogram(x::Vector{Float64}) = collect(0:length(x)-1)' * x / sum(x)

"""
    mean_histogram(x::Vector{Array})

Calculates the mean of multiple histograms.

# Arguments
- `x`: A vector of histograms, where each histogram is represented as a vector of Float64 values.

# Description
This function calculates the mean of each histogram in the vector `x` by calling the `mean_histogram` function for each histogram.

# Returns
- `Vector{Float64}`: A vector of means for each histogram in the input vector.
"""
function mean_histogram(x::Vector{Array})
    y = Vector{Float64}(undef, length(x))
    for i in eachindex(x)
        y[i] = mean_histogram(x[i])
    end
    return y
end

"""
    Difference_Zscore(x1, x2, sig1, sig2)

Calculates the Z-score for the difference between two sets of values.

# Arguments
- `x1`: The first set of values.
- `x2`: The second set of values.
- `sig1`: The standard deviations corresponding to the first set of values.
- `sig2`: The standard deviations corresponding to the second set of values.

# Description
This function calculates the Z-score for the difference between two sets of values `x1` and `x2`. The Z-score is computed as the difference between the values divided by the square root of the sum of the variances.

# Returns
- `Float64`: The Z-score for the difference between the two sets of values.
"""
Difference_Zscore(x1, x2, sig1, sig2) = (x1 - x2) / sqrt(sig1^2 + sig2^2)

"""
    m2_histogram(x)

Calculates the second moment of a histogram.

# Arguments
- `x`: A vector representing the histogram.

# Description
This function calculates the second moment of a histogram `x` by taking the dot product of the squared bin indices and the histogram values, then dividing by the sum of the histogram values.

# Returns
- `Float64`: The second moment of the histogram.
"""
m2_histogram(x) = (collect(0:length(x)-1) .^ 2)' * x / sum(x)

"""
    var_histogram(x)

Calculates the variance of a histogram.

# Arguments
- `x`: A vector representing the histogram.

# Description
This function calculates the variance of a histogram `x` by subtracting the square of the mean of the histogram from the second moment of the histogram.

# Returns
- `Float64`: The variance of the histogram.
"""
function var_histogram(x)
    m2_histogram(x) - mean_histogram(x)^2
end

"""
    moment_histogram(x, n)

Calculates the nth central moment of a histogram.

# Arguments
- `x`: A vector representing the histogram.
- `n`: The order of the moment to calculate.

# Description
This function calculates the nth central moment of a histogram `x` by taking the dot product of the nth power of the deviations from the mean and the histogram values, then dividing by the sum of the histogram values.

# Returns
- `Float64`: The nth central moment of the histogram.
"""
function moment_histogram(x, n)
    y = collect(0:length(x)-1)
    v = (y .- mean_histogram(x)) .^ n
    v' * x / sum(x)
end

"""
    factorial_moment(h::Vector, n)

Calculates the nth factorial moment of a histogram.

# Arguments
- `h`: A vector representing the histogram.
- `n`: The order of the factorial moment to calculate.

# Description
This function calculates the nth factorial moment of a histogram `h` by iterating through the histogram values and computing the product of the bin index and the histogram value for each bin.

# Returns
- `Float64`: The nth factorial moment of the histogram.
"""
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

"""
    moment_param_estimates(h)

Estimates parameters based on the moments of a histogram.

# Arguments
- `h`: A vector representing the histogram.

# Description
This function calculates the first three factorial moments of a histogram `h` and uses them to estimate parameters `lambda`, `mu`, and `nu`.

# Returns
- `Tuple{Float64, Float64, Float64}`: The estimated parameters `lambda`, `mu`, and `nu`.
"""
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
    tstat_2sample(x1, x2)

Computes the t-statistic for the difference between the means of two histograms.

# Arguments
- `x1`: The first histogram.
- `x2`: The second histogram.

# Description
This function computes the t-statistic for the difference between the means of two histograms `x1` and `x2`. The t-statistic is calculated as the difference between the means divided by the square root of the sum of the variances divided by their respective lengths.

# Returns
- `Float64`: The t-statistic for the difference between the means of the two histograms.
"""
tstat_2sample(x1, x2) = (mean_histogram(x1) - mean_histogram(x2)) / (sqrt(var_histogram(x1) / length(x1) + var_histogram(x2) / length(x2)))


"""
    log_2sample(x1, x2)

Computes the log of the ratio of the means of two histograms.

# Arguments
- `x1`: The first histogram.
- `x2`: The second histogram.

# Description
This function computes the log of the ratio of the means of two histograms `x1` and `x2`.

# Returns
- `Float64`: The log of the ratio of the means of the two histograms.
"""
log_2sample(x1, x2) = log(mean_histogram(x1)) - log(mean_histogram(x2))

"""
    delta_2sample(x1, x2)

Computes the difference of the means of two histograms.

# Arguments
- `x1`: The first histogram.
- `x2`: The second histogram.

# Description
This function computes the difference of the means of two histograms `x1` and `x2`.

# Returns
- `Float64`: The difference of the means of the two histograms.
"""
delta_2sample(x1, x2) = mean_histogram(x1) - mean_histogram(x2)

"""
    delta_2frac(x1, x2)

Computes the normalized difference of the means of two histograms.

# Arguments
- `x1`: The first histogram.
- `x2`: The second histogram.

# Description
This function computes the difference of the means of two histograms `x1` and `x2`, normalized by the mean of the first histogram `x1`.

# Returns
- `Float64`: The normalized difference of the means of the two histograms.
"""
delta_2frac(x1, x2) = delta_2sample(x1, x2) / mean_histogram(x1)

"""
    mediansmooth(xin, window)

Applies median smoothing to a histogram over a specified window.

# Arguments
- `xin`: The input histogram.
- `window`: The size of the smoothing window.

# Description
This function applies median smoothing to the input histogram `xin` over the specified window size. It replaces each value in the histogram with the median of the values in the window centered at that value.

# Returns
- `Vector{Float64}`: The smoothed histogram.
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
    residenceprob_G(r::Vector, G::Int)

Calculates the residence probability of G states given the rates.

# Arguments
- `r`: A vector representing the rates.
- `G`: The number of states.

# Description
This function calculates the residence probability of G states given the rates. It initializes the residence probability array `Gss` and iteratively computes the residence probabilities for each state using the rates `r`.

# Returns
- `Array{Float64,2}`: A 2D array containing the residence probabilities for each state.
"""
function residenceprob_G(r::Vector, G::Int)
    n = G - 1
    Gss = Array{Float64,2}(undef, 1, n + 1)
    Gss[1, 1] = 1.0
    for k in 1:n
        Gss[1, k+1] = Gss[1, k] * r[2*k-1] / r[2*k]
    end
    Gss ./= sum(Gss)
end

"""
    residenceprob_G_dataframe(r::Vector, G::Int)

Creates a DataFrame with residence probabilities for G states for a single gene.

# Arguments
- `r`: Vector of rates
- `G`: Number of G states

# Returns
- DataFrame with columns for each G state and their residence probabilities
"""
function residenceprob_G_dataframe(r::Vector, G::Int)
    # Calculate residence probabilities
    probs = vec(residenceprob_G(r, G))

    # Create DataFrame with state labels as columns
    df = DataFrame()
    for j in 1:G
        df[!, "G$j"] = [probs[j]]
    end

    return df
end

function write_joint_residence_prob(outfile, datapath, datacond, interval::Float64, r::Vector, transitions, G, R, S, insertstep, start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)
    traces = read_tracefiles(datapath, datacond, start, stop)
    # r = readrates(file, get_row("median"))
    p0 = make_p0_coupled(traces, interval, r, transitions, G, R, S, insertstep, probfn, noiseparams, splicetype, state, hierarchical, coupling, grid, zeromedian)
    Gjoint, Rjoint = joint_residence_prob(p0, G, R, S, insertstep, coupling)
    df = joint_residence_prob_dataframe(Gjoint, Rjoint)
    CSV.write(outfile, df)
end

function make_p0_coupled(traces, interval, rin, transitions, G, R, S, insertstep, probfn=prob_Gaussian, noiseparams=4, splicetype="", state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)
    _, model = make_trace_datamodel(traces, interval, rin, transitions, G, R, S, insertstep, probfn, noiseparams, splicetype, state, hierarchical, coupling, grid, zeromedian)
    rates, noiseparams, couplingStrength = prepare_rates(get_param(model), model)
    Qtr = make_mat_TC(model.components, rates, couplingStrength)
    p0 = normalized_nullspace(Qtr)
    return p0
end

function joint_residence_prob(p0, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, coupling::Tuple)
     Gjoint = zeros(G...)
     Rjoint = zeros(2,2)
     states = []
    for i in eachindex(p0)
        # Decode the state for each unit
        units = unit_state(i, G, R, S, coupling[1])
        # states is a tuple of (g, z, zdigits, r) for each unit
        states = inverse_state(units, G, R, S, insertstep, any)
        g_indices = [s[1] for s in states]  # Extract G state for each unit
        r_indices = [s[4] for s in states] .+ 1
        Gjoint[g_indices...] += p0[i]
        Rjoint[r_indices...] += p0[i]
    end
    # return states
    return Gjoint ./ sum(Gjoint), Rjoint ./ sum(Rjoint)
end

"""
    joint_residence_prob_dataframe(Gjoint, Rjoint)

Creates a DataFrame with joint residence probabilities for G states and R states.

# Arguments
- `Gjoint`: N-dimensional array of joint G state probabilities
- `Rjoint`: 2D array of joint R state probabilities

# Returns
- DataFrame with columns for each G state, R state, and their joint probabilities
"""
function joint_residence_prob_dataframe(Gjoint, Rjoint)
    # Initialize DataFrame
    df = DataFrame()
    
    # Add G state columns
    G_dims = size(Gjoint)
    n_genes = length(G_dims)
    
    # Create all possible combinations of G states
    G_indices = CartesianIndices(Gjoint)
    G_combinations = [Tuple(idx) for idx in G_indices]
    
    # Add G state columns
    for i in 1:n_genes
        df[!, "G$(i)"] = [comb[i] for comb in G_combinations]
    end
    
    # Add G joint probability
    df[!, :G_Probability] = [Gjoint[idx...] for idx in G_indices]
    
    # Add R state columns
    R_dims = size(Rjoint)
    R_indices = CartesianIndices(Rjoint)
    R_combinations = [Tuple(idx) for idx in R_indices]
    
    # Add R state columns with OFF/ON labels
    for i in 1:length(R_dims)
        df[!, "Reporter$(i)"] = [comb[i] == 1 ? "OFF" : "ON" for comb in R_combinations]
    end
    
    # Add R joint probability
    df[!, :Reporter_Probability] = [Rjoint[idx...] for idx in R_indices]
    
    return df
end

"""
    residenceprob_G_dataframe(r::Vector, G::Tuple, nrates::Vector{Int})

Creates a DataFrame with residence probabilities for G states across multiple genes.

# Arguments
- `r`: Vector of rates
- `G`: Tuple containing number of G states for each gene
- `nrates`: Vector containing number of rates for each gene

# Returns
- DataFrame with columns for each G state and their residence probabilities
"""
function residenceprob_G_dataframe(r::Vector, G::Tuple, nrates::Vector{Int})
    # Initialize DataFrame
    df = DataFrame()
    k = 1

    # Calculate residence probabilities for each gene using the Int version
    for (i, g) in enumerate(G)
        # Get single gene probabilities using the Int version
        single_df = residenceprob_G_dataframe(r[k:k+nrates[i]-1], g)

        # Rename columns to include gene number
        for col in names(single_df)
            df[!, "$(col)_$(i)"] = single_df[!, col]
        end

        k += nrates[i]
    end

    return df
end

"""
    splicesiteusage(model::AbstractGRSMmodel)
    splicesiteusage(r::Vector, n::Int, nr::Int)

Calculates the splice site usage for a GRSM model or given rates.

# Arguments
- `model`: An instance of `AbstractGRSMmodel`.
- `r`: Rates vector.
- `G`: Number of states
- `S`: Number of splice sites.

# Methods

## `splicesiteusage(model::AbstractGRSMmodel)`

Calculates the splice site usage for a GRSM model.

## `splicesiteusage(r::Vector, n::Int, nr::Int)`

Calculates the splice site usage given rates, number of states, and number of splice sites.

# Returns
- `Vector{Float64}`: Splice site usage probabilities.
"""
splicesiteusage(model::AbstractGRSMmodel) = splicesiteusage(model.rates, model.G, model.R)

function splicesiteusage(r::Vector, G::Int, S::Int)
    n = G - 1
    nr = S
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
    onstate_prob(r::Vector, components::TComponents, reporter_per_state)
    onstate_prob(param, model::AbstractGeneTransitionModel)
    onstate_prob(model::AbstractGeneTransitionModel)

Calculates the probability of being in the "on" state for a given model or rates.

# Arguments
- `r`: Rates vector.
- `components`: An instance of `TComponents`.
- `reporter_per_state`: A vector indicating the reporter per state.
- `param`: Parameters for the model.
- `model`: An instance of `AbstractGeneTransitionModel`.

# Methods

## `onstate_prob(r::Vector, components::TComponents, reporter_per_state)`

Calculates the probability of being in the "on" state given the rates, components, and reporter per state.

## `onstate_prob(param, model::AbstractGeneTransitionModel)`

Calculates the probability of being in the "on" state given the parameters and model.

## `onstate_prob(model::AbstractGeneTransitionModel)`

Calculates the probability of being in the "on" state for a given model.

# Returns
- `Tuple{Float64, Vector{Float64}}`: The probability of being in the "on" state and the steady-state probabilities.
"""
function onstate_prob(r, components::TComponents, reporter_per_state)
    Qtr = make_mat(components.elementsT, r, components.nT)
    p0 = normalized_nullspace(Qtr)
    return sum(p0[reporter_per_state.>0]), p0
end

onstate_prob(r, model::AbstractGeneTransitionModel) = onstate_prob(r, model.components, model.reporter.per_state)

# onstate_prob(param, model::AbstractGeneTransitionModel) = onstate_prob(get_rates(param, model), model.components, model.reporter.per_state)

onstate_prob(model::AbstractGeneTransitionModel) = onstate_prob(model.rates, model.components, model.reporter.per_state)

"""
Obsolete function. needs to be updated.

    burstoccupancy(model::AbstractGRSMmodel)
    burstoccupancy(n::Int, nr::Int, r::Vector)

Calculates the burst size distribution of a GRS model for total pre-mRNA occupancy and unspliced (visible) occupancy.

# Description
This function calculates the burst size distribution for a given GRS model. It computes the total pre-mRNA occupancy and unspliced (visible) occupancy distributions by marginalizing over the conditional steady-state distribution.

# Arguments
- `model`: An instance of `AbstractGRSMmodel`.
- `n`: The number of states.
- `nr`: The number of splice sites.
- `r`: A vector representing the rates.

# Methods

## `burstoccupancy(model::AbstractGRSMmodel)`

Calculates the burst size distribution for a given GRSM model.

## `burstoccupancy(n::Int, nr::Int, r::Vector)`

Calculates the burst size distribution given the number of states, number of splice sites, and rates.

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: The total pre-mRNA occupancy and unspliced (visible) occupancy distributions.
"""
burstoccupancy(model::AbstractGRSMmodel) = burstoccupancy(model.G - 1, model.R, model.rates)

function burstoccupancy(n::Int, nr::Int, r::Vector)
    T = mat_GSR_T(r, n, nr)
    pss = normalized_nullspace(T)
    Rss = zeros(nr)
    Rssvisible = zeros(nr)
    ssf = zeros(nr)
    asum = 0
    for w = 1:nr
        for i in 1:n+1, z = 1:3^nr
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

"""
    model2_mean(ron, roff, eject, decay, nalleles)
    model2_mean(ron, roff, eject, decay)

Calculates the mean mRNA level for a two-state model.

# Arguments
- `ron`: The rate of switching to the "on" state.
- `roff`: The rate of switching to the "off" state.
- `eject`: The rate of mRNA ejection.
- `decay`: The rate of mRNA decay.
- `nalleles`: The number of alleles (optional).

# Methods

## `model2_mean(ron, roff, eject, decay, nalleles)`

Calculates the mean mRNA level for a two-state model with multiple alleles.

## `model2_mean(ron, roff, eject, decay)`

Calculates the mean mRNA level for a two-state model with a single allele.

# Returns
- `Float64`: The mean mRNA level.
"""
model2_mean(ron, roff, eject, decay, nalleles) = nalleles * model2_mean(ron, roff, eject, decay)

model2_mean(ron, roff, eject, decay) = ron / (ron + roff) * eject / decay

"""
    model2_variance(ron, roff, eject, decay, nalleles)
    model2_variance(ron, roff, eject, decay)

Calculates the variance of mRNA levels for a two-state model.

# Arguments
- `ron`: The rate of switching to the "on" state.
- `roff`: The rate of switching to the "off" state.
- `eject`: The rate of mRNA ejection.
- `decay`: The rate of mRNA decay.
- `nalleles`: The number of alleles (optional).

# Methods

## `model2_variance(ron, roff, eject, decay, nalleles)`

Calculates the variance of mRNA levels for a two-state model with multiple alleles.

## `model2_variance(ron, roff, eject, decay)`

Calculates the variance of mRNA levels for a two-state model with a single allele.

# Returns
- `Float64`: The variance of mRNA levels.
"""
model2_variance(ron, roff, eject, decay, nalleles) = nalleles * model2_variance(ron, roff, eject, decay)

model2_variance(ron, roff, eject, decay) = ron / (ron + roff) * eject / decay + ron * roff / (ron + roff)^2 * eject^2 / decay / (ron + roff + decay)


"""
    multi_tau_covariance(data1::Vector{Float64}, data2::Vector{Float64}=Float64[]; 
                         m::Int=16, max_lag::Int=length(data1))

Compute auto/cross-covariance or correlation using multi-tau algorithm.

# Arguments
- `data1::Vector{Float64}`: First time series
- `data2::Vector{Float64}`: Optional second time series for cross-correlation
- `m::Int=16`: Points per level
- `max_lag::Int`: Maximum lag time (default: length of data)
- `normalize::Bool=false`: Normalize by variance to obtain correlation values (default: false)

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: (lag times, correlation values)
"""
function multi_tau_covariance(data1::Vector{Float64}, data2::Vector{Float64}=Float64[];
    m::Int=16, max_lag::Int=length(data1), normalize=false)

    n = length(data1)
    is_auto = isempty(data2)
    data2 = is_auto ? data1 : data2

    # Determine number of levels needed
    max_level = floor(Int, log2(max_lag / m)) + 1
    total_points = m * (2^max_level - 1)
    num_lags = min(total_points, max_lag)

    # Initialize output arrays
    lags = zeros(Float64, num_lags)
    corr = zeros(Float64, num_lags)
    counts = zeros(Int, num_lags)

    # Calculate means for normalization
    mean1 = mean(data1)
    mean2 = mean(data2)

    # First level correlations
    for τ in 1:min(m, num_lags)
        lags[τ] = τ
        for i in 1:(n-τ)
            corr[τ] += (data1[i] - mean1) * (data2[i+τ] - mean2)
            counts[τ] += 1
        end
    end

    # Higher levels with increasing spacing
    lag_idx = m + 1
    for level in 1:max_level-1
        spacing = 2^level

        for τ in (m*spacing):(2*m*spacing-1)
            if τ > num_lags
                break
            end

            lags[lag_idx] = τ
            for i in 1:(n-τ)
                corr[lag_idx] += (data1[i] - mean1) * (data2[i+τ] - mean2)
                counts[lag_idx] += 1
            end
            lag_idx += 1
        end
    end
    if normalize
        # Normalize correlations
        var1 = var(data1)
        var2 = is_auto ? var1 : var(data2)
        norm_factor = sqrt(var1 * var2)
    else
        norm_factor = 1.0
    end

    for i in 1:num_lags
        if counts[i] > 0
            corr[i] /= (counts[i] * norm_factor)
        end
    end


    # Remove unused lag times
    valid_idx = counts .> 0
    return lags[valid_idx], corr[valid_idx]
end

"""
    make_histogram(r)

Creates a histogram from a vector of data.

# Arguments
- `r`: A vector of data values.

# Description
This function creates a histogram from the provided vector of data values `r`. It counts the occurrences of each unique value in the vector and returns a vector representing the histogram.

# Returns
- `Vector{Int}`: A vector representing the histogram, where each element corresponds to the count of a unique value in the input vector.
"""
function make_histogram(r; normalize=false)
    c = Int.(r)
    nhist = maximum(c)
    h = zeros(Int, nhist + 1)
    for i in c
        h[i+1] += 1
    end
    if normalize
        return h / sum(h)
    else
        return h
    end
end

"""
    make_histograms(folder::String, file::String, label::String)

Creates histograms from cell counts data in a file and saves them to a specified folder.

# Arguments
- `folder`: The folder where the histograms will be saved.
- `file`: The file containing the data.
- `label`: A label to be appended to the histogram filenames.

# Description
This function reads data from the specified file, creates histograms for each row of data, and saves the histograms to the specified folder. Each histogram is saved with a filename that includes the row identifier and the provided label.

# Returns
- `Nothing`: The function performs the histogram creation and saving but does not return a value.
"""


function make_rna_histograms(folder, count_matrix::Matrix, label)
    if ~ispath(folder)
        mkpath(folder)
    end
    for r in eachrow(count_matrix)
        f = open("$folder/$(r[1])_$label.txt", "w")
        a = r[2:end]
        a = a[(a.!="").&(a.!="NA")]
        h = make_histogram(Int.(a) .+ 1)
        writedlm(f, h)
        close(f)
    end
end

function make_rna_histograms(folder, file::String, label::String)
    count_matrix, _ = read_cell_counts(file)   
    make_rna_histograms(folder, count_matrix, label)
end



function get_cell_counts(count_matrix::Matrix)
    if typeof(count_matrix[1, 1]) <: AbstractString
        c = sum(count_matrix[:, 2:end], dims=1)
    else
        c = sum(count_matrix, dims=1)
    end
    S = sum(c)
    return vec(c) / maximum(c), S
end

"""
    make_gene_counts(folder, gene::String, counts::Vector, yieldfactor::Vector, label::String)

Creates a file containing gene counts and yield factors for a single gene.

# Arguments
- `folder::String`: Output directory where the file will be saved
- `gene::String`: Name of the gene
- `counts::Vector`: Vector of RNA counts across cells
- `yieldfactor::Vector`: Vector of yield factors for normalization
- `label::String`: Label to append to the output filename

# Description
Creates a file named `{gene}_{label}.txt` containing two columns:
1. RNA counts for each cell
2. Corresponding yield factors for normalization

# Example
```julia
make_gene_counts("output", "GENE1", [1,2,3,4], [0.5,0.5,0.5,0.5], "control")
```
"""
function make_gene_counts(folder, gene::String, counts::Vector, yieldfactor::Vector, label::String)
    f = open("$folder/$(gene)_$label.txt", "w")
    writedlm(f, [counts yieldfactor])
    close(f)
end

"""
    make_gene_counts(folder, count_matrix::Matrix, yieldfactor::Vector, label::String)

Processes a count matrix to create individual files for each gene containing counts and yield factors.

# Arguments
- `folder::String`: Output directory where files will be saved
- `count_matrix::Matrix`: Matrix where each row represents a gene and columns represent cells
- `yieldfactor::Vector`: Vector of yield factors for normalization
- `label::String`: Label to append to output filenames

# Description
For each gene in the count matrix:
1. Creates a file named `{gene}_{label}.txt`
2. Writes two columns:
   - RNA counts for each cell (filtered to remove empty and NA values)
   - Corresponding yield factors for normalization

# Example
```julia
counts = [1 2 3; 4 5 6]  # 2 genes, 3 cells
yield = [0.5, 0.5, 0.5]
make_gene_counts("output", counts, yield, "control")
```
"""
function make_gene_counts(folder, count_matrix::Matrix, yieldfactor::Vector, label::String)
    if ~ispath(folder)
        mkpath(folder)
    end
    for r in eachrow(count_matrix)
        gene = replace(replace(string(r[1]), "/" => "-"), "\\" => "-")
        counts = r[2:end]
        counts = counts[(counts.!="").&(counts.!="NA")]
        make_gene_counts(folder, gene, counts, yieldfactor, label)
    end
end

"""
    make_gene_counts(folder, file, label)

Reads a count matrix from a file and creates individual files for each gene.

# Arguments
- `folder::String`: Output directory where files will be saved
- `file::String`: Path to input file containing count matrix
- `label::String`: Label to append to output filenames

# Description
1. Reads count matrix from input file
2. Calculates yield factors from cell counts
3. Creates individual files for each gene containing:
   - RNA counts for each cell
   - Corresponding yield factors

# Example
```julia
make_gene_counts("output", "counts.csv", "control")
```
"""
function make_gene_counts(folder, file, label)
    count_matrix, _ = readdlm(file, header=true)
    yieldfactor, _ = get_cell_counts(count_matrix)
    make_gene_counts(folder, count_matrix, yieldfactor, label)
end





"""
    expv(v::Array)

Calculates the element-wise exponential of an array.

# Arguments
- `v`: An array of numerical values.

# Description
This function calculates the element-wise exponential of the provided array `v`. It returns a new array where each element is the exponential of the corresponding element in `v`.

# Returns
- `Array`: An array containing the element-wise exponential of the input array.
"""
expv(v::Array) = exp.(v)

"""
    logv(v::Array)

Calculates the element-wise natural logarithm of an array.

# Arguments
- `v`: An array of numerical values.

# Description
This function calculates the element-wise natural logarithm of the provided array `v`. It returns a new array where each element is the natural logarithm of the corresponding element in `v`.

# Returns
- `Array`: An array containing the element-wise natural logarithm of the input array.
"""
logv(v::Array) = log.(v)


"""
    log_shift(v::Float64, a::Float64)
    log_shift(v::Array, a::Float64)

Calculates the natural logarithm of a value or array after applying a shift.

# Arguments
- `v`: A numerical value or an array of numerical values.
- `a`: The shift value to be added before taking the logarithm.

# Methods

## `log_shift(v::Float64, a::Float64)`

Calculates the natural logarithm of a single value after applying a shift.

## `log_shift(v::Array, a::Float64)`

Calculates the element-wise natural logarithm of an array after applying a shift.

# Returns
- `Float64` or `Array`: The natural logarithm of the shifted value or array.
"""
log_shift(v::Float64, a::Float64) = log(v + a)

log_shift(v::Array, a::Float64) = log.(v .+ a)

"""
    log_shift1(v)

Calculates the natural logarithm of a value or array after applying a shift of 1.0.

# Arguments
- `v`: A numerical value or an array of numerical values.

# Description
This function calculates the natural logarithm of the provided value or array `v` after applying a shift of 1.0. It uses the `log_shift` function with a shift value of 1.0.

# Returns
- `Float64` or `Array`: The natural logarithm of the shifted value or array.
"""
log_shift1(v) = log_shift(v, 1.0)

"""
    invlog_shift1(v::Float64)
    invlog_shift1(v::Array)

Calculates the inverse of the natural logarithm after applying a shift of 1.0.

# Arguments
- `v`: A numerical value or an array of numerical values.

# Methods

## `invlog_shift1(v::Float64)`

Calculates the inverse of the natural logarithm for a single value after applying a shift of 1.0.

## `invlog_shift1(v::Array)`

Calculates the element-wise inverse of the natural logarithm for an array after applying a shift of 1.0.

# Returns
- `Float64` or `Array`: The inverse of the natural logarithm of the shifted value or array.
"""
invlog_shift1(v::Float64) = exp(v) - 1
invlog_shift1(v::Array) = exp.(v) .- 1

"""
    logit(x::Float64)
    logit(x::Array)

Applies the logit transform to a value or array.

# Arguments
- `x`: A numerical value or an array of numerical values.

# Methods

## `logit(x::Float64)`

Applies the logit transform to a single value.

## `logit(x::Array)`

Applies the logit transform element-wise to an array.

# Description
The logit transform is defined as `log(x) - log(1 - x)`.

# Returns
- `Float64` or `Array`: The logit-transformed value or array.
"""
logit(x::Float64) = log(x) - log(1 - x)
logit(x::Array) = log.(x) .- log.(1 .- x)

"""
    invlogit(x::Float64)
    invlogit(x::Array)

Applies the inverse logit transform to a value or array.

# Arguments
- `x`: A numerical value or an array of numerical values.

# Methods

## `invlogit(x::Float64)`

Applies the inverse logit transform to a single value.

## `invlogit(x::Array)`

Applies the inverse logit transform element-wise to an array.

# Description
The inverse logit transform is defined as `1 / (1 + exp(-x))`.

# Returns
- `Float64` or `Array`: The inverse logit-transformed value or array.
"""
invlogit(x::Float64) = 1 / (1 + exp(-x))
invlogit(x::Array) = 1 ./ (1 .+ exp.(-x))

"""
    logsumexp(u::Array, v::Array)
    logsumexp(u::Float64, v::Float64)
    logsumexp(v::Array)

Returns the log of the sum of exponentials of the inputs.

# Arguments
- `u`: An array or a single floating-point value.
- `v`: An array or a single floating-point value.

# Methods

## `logsumexp(u::Array, v::Array)`

Calculates the log of the sum of exponentials for two arrays element-wise.

## `logsumexp(u::Float64, v::Float64)`

Calculates the log of the sum of exponentials for two floating-point values.

## `logsumexp(v::Array)`

Calculates the log of the sum of exponentials for an array.

# Returns
- `Array` or `Float64`: The log of the sum of exponentials of the inputs.
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

function logsumexp(v::Array)
    w = maximum(v)
    if w == -Inf
        return -Inf
    else
        return w .+ log(sum(exp.(v .- w)))
    end
end
"""
    kronecker_delta(i, j)

Computes the Kronecker delta of two integers.

# Arguments
- `i`: The first integer.
- `j`: The second integer.

# Returns
- `Int`: Returns 1 if `i` equals `j`, otherwise returns 0.
"""
kronecker_delta(i, j) = i == j ? 1 : 0

"""
    meantime(r::Vector)

Calculates the mean time given a vector of rates.

# Arguments
- `r`: A vector of rates.

# Description
This function calculates the mean time by taking the reciprocal of each rate in the vector `r` and summing them up.

# Returns
- `Float64`: The mean time.
"""
function meantime(r::Vector)
    sum(1 ./ r)
end

"""
    mean_elongationtime(r, model)

Calculates the mean elongation time for a given model and rates.

# Arguments
- `r`: The rates.
- `model`: The model, which includes the elongation rates and states.

# Description
This function calculates the mean elongation time for a given model and rates. It uses the elongation rates and states specified in the model to compute the mean time required for elongation.

# Returns
- `Float64`: The mean elongation time.
"""
function mean_elongationtime(r, transitions, R::Int)
    if R > 0
        n = length(transitions)
        return meantime(r[n+2:n+R+1])
    else
        return 1.0
    end
end

function mean_elongationtime(r, transitions, R::Tuple)
    m = Vector{Float64}(undef, length(R))
    for i in eachindex(R)
        m[i] = mean_elongationtime(r, transitions, R[i])
    end
    return m
end
"""
    make_array(v)

Converts a container of arrays into a single array.

# Arguments
- `v`: A container of arrays.

# Description
This function concatenates the arrays in the container `v` into a single array.

# Returns
- `Array`: A single array containing the elements of the input container.

"""
function make_array(v)
    vconcat = Float64[]
    for v in v
        vconcat = [vconcat; v]
    end
    return vconcat
end

"""
    makestring(v)

convert a vector of strings into a single string

# Arguments
- `v`: A vector of strings.

# Description
This function concatenates the strings in the vector `v` into a single string.

# Returns
- `String`: A single string containing the elements of the input vector.

"""
function makestring(v)
    s = ""
    for i in v
        s *= i
    end
    return s
end

function grid_distance(i, j, Nx, Ny)
    if Nx == 0 || Ny == 0
        return 1.0
    else
        xi = rem(i - 1, Nx)
        yi = div(i - 1, Ny)
        xj = rem(j - 1, Nx)
        yj = div(j - 1, Ny)
        return sqrt((xi - xj)^2 + (yi - yj)^2)
    end
end

"""
    grid_distance(i, j, Ngrid)

Calculate the Euclidean distance between two points in a grid.

# Arguments
- `i`: The index of the first point.
- `j`: The index of the second point.
- `Ngrid`: The number of points along one dimension of the grid.

# Description
This function calculates the Euclidean distance between two points `i` and `j` in a grid with `Ngrid` points along one dimension. The points are assumed to be arranged in a 2D grid.

# Returns
- `Float64`: The Euclidean distance between the two points.

"""
function grid_distance(i, j, Ngrid)
    grid_distance(i, j, Ngrid, Ngrid)
end
# function grid_distance(i, j, Ngrid)
#     if Ngrid == 0
#         return 1.0
#     else
#         xi = rem(i - 1, Ngrid)
#         yi = div(i - 1, Ngrid)
#         xj = rem(j - 1, Ngrid)
#         yj = div(j - 1, Ngrid)
#         return sqrt((xi - xj)^2 + (yi - yj)^2)
#     end
# end

"""
    find_top_two(x::AbstractArray{T}) where T <: Number

Finds the largest and second largest values in an array.

# Arguments
- `x`: An array of numbers.

# Returns
- `Tuple{T, T, Int}`: (largest value, second largest value, index of largest value)

# Throws
- `ArgumentError`: If array has fewer than 2 elements.
"""
function find_top_two(x::AbstractArray{T}) where {T<:Number}
    length(x) < 2 && throw(ArgumentError("Array must have at least 2 elements"))

    largest = second = typemin(T)
    largest_idx = 0

    for (i, val) in enumerate(x)
        if val > largest
            second = largest
            largest = val
            largest_idx = i
        elseif val > second && val < largest
            second = val
        end
    end

    return largest, second, largest_idx
end

"""
    replace_outlier!(x::AbstractArray{T}) where T <: Number

Replaces the largest value with second largest if it's more than double in magnitude.

# Arguments
- `x`: Array to check and potentially modify.

# Returns
- `Bool`: true if a replacement was made, false otherwise.
"""
function replace_outlier!(x::AbstractArray{T}) where {T<:Number}
    largest, second, idx = find_top_two(x)
    if abs(largest) > 2 * abs(second)
        x[idx] = second
        return true
    end
    return false
end



