using Random

rng = Xoshiro(13)

# Hamiltonian MC - no gradients, no symplectic integration, no tuning
function H(f::Function, g::Function, z::Tuple)
    f_ = f(z[1])
    g_ = g(f_)
    log_g = log(g_)
    z_next = z[1] + z[2] * g_
    f_next = f(z_next)
    g_next = g(f_next)
    return f_next + f_ - log_g, -log(g_next) + log_g, (z_next, -z[2] * g_ / g_next)
end

function H_neg(f::Function, g::Function, z::Tuple)
    f_ = f(z[1])
    g_ = g(f_)
    log_g = log(g_)
    z_next = z[1] - z[2] * g_
    f_next = f(z_next)
    g_next = g(f_next)
    return f_next + f_ - log_g, -log(g_next) + log_g, (z_next, -z[2] * g_ / g_next)
end

function Steffensen_IHMC(f::Function, a::Number, beta::Float64;
    p_range::Float64=5.0, burn_in=3000,
    tol::AbstractFloat=1e-10, maxiter::Integer=20000,
    print_iter::Integer=45)
    # Sigmoid-like function for mapping fx
    sigm = fx -> exp(beta * fx) # Must be a monotonic positive function

    # Hamiltonian functions
    Hpq = z -> H(f, sigm, z)
    Hpq_neg = z -> H_neg(f, sigm, z)

    # Initialization
    it_ = 0
    local f1, z, del_H, H_prev
    z = (a, (rand(rng) - 0.5) * p_range)
    H_prev = Hpq(z)

    # Iterative sampling
    while it_ < maxiter
        if it_ % 2 == 1
            H_prev, del_H, z_next = Hpq(z)
        else
            H_prev, del_H, z_next = Hpq_neg(z)
        end

        # Metropolis acceptance step
        accept = 0
        if del_H < -log(rand(rng))
            accept = 1
        end

        # Logging
        if it_ % (maxiter // print_iter) == 0 && it_ > burn_in
            println(it_, " q ", z[1], " p ", z[2], " H ", H_prev,
                " del_H ", del_H, " del_f ", f(z_next[1]) - f(z[1]),
                " accept ", accept)
        end

        # Update state
        if accept == 1
            z = (z_next[1], (rand(rng) - 0.5) * p_range)
        else
            z = (z[1], (rand(rng) - 0.5) * p_range)
        end

        it_ += 1
    end

    return H_prev + del_H, H_prev
end
