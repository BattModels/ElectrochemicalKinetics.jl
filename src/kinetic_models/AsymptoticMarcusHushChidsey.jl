"""
    AsymptoticMarcusHushChidsey(A, λ)
    AsymptoticMarcusHushChidsey(λ)

Computes asymptotic solution to MHC model, as described in Zeng et al.: 10.1016/j.jelechem.2014.09.038d, with a corrction prefactor of kT since there is an error in the nondimensionalization in that work.

If initialized with one argument, assumes this to be the reorganization energy λ and sets the prefactor to 1.0.
"""
struct AsymptoticMarcusHushChidsey{T} <: NonIntegralModel{T}
    A::T
    λ::T
    function AsymptoticMarcusHushChidsey(A, λ)
        ps = consistent_params(Float32.(A), Float32.(λ))
        new{typeof(ps[1])}(ps...)
    end
end

# default prefactor is 1
AsymptoticMarcusHushChidsey(λ) = AsymptoticMarcusHushChidsey(1.0, λ)

function rate_constant(V_app, amhc::AsymptoticMarcusHushChidsey, ox::Bool; T = 298)
    η = (2 * ox - 1) .* V_app ./ (kB * T)
    λ_nondim = amhc.λ / (kB * T)
    a = 1 .+ sqrt.(λ_nondim)
    arg = (λ_nondim .- sqrt.(a .+ η.^2)) ./ (2 .* sqrt.(λ_nondim))
    pref = sqrt.(π .* λ_nondim) ./ (1 .+ exp.(-η))
    return kB * T .* amhc.A .* pref .* erfc.(arg)
end

rate_constant(V_app, amhc::AsymptoticMarcusHushChidsey, ::Val{true}; T=298) = rate_constant(V_app, amhc, true; T=T)
rate_constant(V_app, amhc::AsymptoticMarcusHushChidsey, ::Val{false}; T=298) = rate_constant(V_app, amhc, false; T=T)