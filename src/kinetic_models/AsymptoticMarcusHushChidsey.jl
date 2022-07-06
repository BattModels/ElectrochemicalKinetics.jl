"""
    AsymptoticMarcusHushChidsey(A, λ)
    AsymptoticMarcusHushChidsey(λ)

Computes asymptotic solution to MHC model, as described in Zeng et al.: 10.1016/j.jelechem.2014.09.038d, with a corrction prefactor of kT since there is an error in the nondimensionalization in that work.

If initialized with one argument, assumes this to be the reorganization energy λ and sets the prefactor to 1.0.
"""
struct AsymptoticMarcusHushChidsey{T,S} <: NonIntegralModel
    A::T
    λ::S
end

# default prefactor is 1
AsymptoticMarcusHushChidsey(λ) = AsymptoticMarcusHushChidsey(1.0, λ)

function amhc_f((V, ps), ::Val{ox}; kT = 0.026) where ox
    η = (2 * ox - 1) .* V ./ kT
    λ_nondim = ps[2] / kT
    a = 1 .+ sqrt.(λ_nondim)
    arg = (λ_nondim .- sqrt.(a .+ η.^2)) ./ (2 .* sqrt.(λ_nondim))
    pref = sqrt.(π .* λ_nondim) ./ (1 .+ exp.(-η))
    kT .* ps[1] .* pref .* erfc.(arg)
end

# direct dispatch for net rates
function amhc_f((V, ps); kT = 0.026)
    η = V/ kT
    λ_nondim = ps[2] / kT
    a = 1 .+ sqrt.(λ_nondim)
    arg = (λ_nondim .- sqrt.(a .+ η.^2)) ./ (2 .* sqrt.(λ_nondim))
    pref = sqrt.(π .* λ_nondim) .* tanh.(η/2)
    return abs.(kT * ps[1] .* pref .* erfc.(arg))
end

rate_f(::AsymptoticMarcusHushChidsey) = amhc_f

# and we have to directly dispatch this case too for it to use the custom net rate case...not the most elegant
function rate_constant(V_app, amhc::AsymptoticMarcusHushChidsey; kT = 0.026)
    res = amhc_f.(Iterators.product(V_app, iterate_props(amhc)); kT=kT)
    if size(res) == (1,)
        return res[1]
    else
        return res
    end
end