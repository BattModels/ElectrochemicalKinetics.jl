"""
    AsymptoticMarcusHushChidsey(A, λ)
    AsymptoticMarcusHushChidsey(λ)

Computes asymptotic solution to MHC model, as described in Zeng et al.: 10.1016/j.jelechem.2014.09.038d, with a corrction prefactor of kT since there is an error in the nondimensionalization in that work.

If initialized with one argument, assumes this to be the reorganization energy λ and sets the prefactor to 1.0.
"""
struct AsymptoticMarcusHushChidsey <: KineticModel
    A::Float64
    λ::Float64
end

# default prefactor is 1
AsymptoticMarcusHushChidsey(λ) = AsymptoticMarcusHushChidsey(1.0, λ)

function (amhc::AsymptoticMarcusHushChidsey)(V_app, ox::Bool; kT::Real = 0.026)
    η = (2 * ox - 1) .* V_app ./ kT
    λ_nondim = amhc.λ / kT
    a = 1 + sqrt(λ_nondim)
    arg = (λ_nondim .- sqrt.(a .+ η.^2)) ./ (2 * sqrt(λ_nondim))
    pref = sqrt(π * λ_nondim) ./ (1 .+ exp.(-η))
    return kT * amhc.A .* pref .* erfc.(arg)
end


function (amhc::AsymptoticMarcusHushChidsey)(V_app; kT::Real = 0.026)
    η = V_app ./ kT
    λ_nondim = amhc.λ / kT
    a = 1 + sqrt(λ_nondim)
    arg = (λ_nondim .- sqrt.(a .+ η.^2)) ./ (2 * sqrt(λ_nondim))
    pref = sqrt(π * λ_nondim) .* tanh.(η/2)
    return abs.(kT * amhc.A .* pref .* erfc.(arg))
end