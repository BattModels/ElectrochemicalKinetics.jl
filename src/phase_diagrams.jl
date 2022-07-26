using ElectrochemicalKinetics
using Plots
using NLsolve
using FastGaussQuadrature

# define some constants and default parameters
const room_T = 298
const muoA_default = 0.02
const muoB_default = 0.03
const Ω_default = 0.1

# our familiar thermodynamic functions
h(x;Ω=Ω_default) = @. x*(1-x)*Ω # enthalpy of mixing
s(x) = @. -kB*(x*log(x+eps(Float32)) + (1-x)*log(1-x+eps(Float32))) # entropy per particle...added epsilons in the hopes that things will be less obnoxious at the edges
g_thermo(x; Ω=Ω_default, muoA=muoA_default, muoB=muoB_default, T=room_T) = @. h(x;Ω) - T*s(x)+ muoA*(1-x) + muoB*x # Gibbs free energy per particle
μ_thermo(x; Ω=Ω_default, muoA=muoA_default, muoB=muoB_default, T=room_T) = @. (1-2*x)*Ω + kB*T*log(x/(1-x)) + muoB-muoA # chemical potential

prefactor(x, intercalate::Bool) = intercalate ? (1 .- x) : x

"""
µ and g with kinetic constributions, can be modeled using any <:KineticModel object

These functions return single-argument functions (to easily use common-tangent function below while
still being able to swap out model parameters by calling "function-builders" with different arguments).
"""
function µ_kinetic(I, km::KineticModel; intercalate=true, kwargs...)
    thermo_term(x) = μ_thermo(x; kwargs...)
    μ(x::Real) = thermo_term(x) .+ overpotential(I, prefactor(x, intercalate)*km)
    μ(x::AbstractVector) = thermo_term(x) .+ overpotential(I, prefactor(x, intercalate).*Ref(km))
    return μ
end

function g_kinetic(I, km::KineticModel; intercalate=true, kwargs...)
    thermo_term(x) = g_thermo(x; kwargs...)
    #TODO: gradient of this term is just value of overpotential(x)
    function kinetic_term(x)
        f(x) = ElectrochemicalKinetics.overpotential(I, prefactor(x, intercalate) * km)
        n, w = ElectrochemicalKinetics.scale_coarse(zero.(x), x)
        map((w, n) -> sum(w .* f(n)), eachcol(w), eachcol(n))
    end
    g(x) = thermo_term(x) .+ kinetic_term(vec(x))
    g(x::Real) = thermo_term(x) + kinetic_term(x)[1]
    return g
end

# zeros of this function correspond to pairs of x's satisfying the common tangent condition for a given µ function
# case where we just feed in two points (x should be of length two)
function common_tangent(x::Vector, I, km::KineticModel; intercalate=true, kwargs...)
    g = g_kinetic(I, km; intercalate=intercalate, kwargs...)
    µ = µ_kinetic(I, km; intercalate=intercalate, kwargs...)
    [(g(x[2]) - g(x[1]))/(x[2] - x[1]) .- μ(x[1]), μ(x[2])-μ(x[1])]
end

# TODO: see if we can speed this up with gradients? And/or if it's even needed for integral cases
function find_phase_boundaries(I, km::KineticModel; guess=[0.05, 0.95], intercalate=true, verbose=false, kwargs...)
    
    function myct!(storage, x)
        res = common_tangent(x, I, km; intercalate=intercalate, kwargs...)
        storage[1] = res[1]
        storage[2] = res[2]
    end
    x1 = nlsolve(myct!, guess, show_trace=verbose, ftol=1e-6, xtol=1e-6)
    x1.zero
end

"""
    phase_diagram(km; I_step=1, I_max=Inf, verbose=false, intercalate=true, kwargs...)

Construct electrochemical phase diagram for the model `km` as a function of composition and current, in steps of `I_step`.
    
Note that appropriate values of `I_step` depend strongly on the prefactor of your model. For example, for `ButlerVolmer` with a prefactor of 1, 
"""
function phase_diagram(km::KineticModel; I_step=1, I_max=Inf, verbose=false, intercalate=true, kwargs...)
    I = 0
    pbs_here = find_phase_boundaries(I, km; intercalate=intercalate, kwargs...)
    pbs = pbs_here'
    I_vals = [0]
    while abs(pbs_here[2] - pbs_here[1]) > 1e-3 && I < I_max
        I = I + I_step
        if verbose
            println("Solving at I=", I, "...")
        end
        try
            pbs_here = find_phase_boundaries(I, km; intercalate=intercalate, guess=pbs_here, kwargs...)
            pbs = vcat(pbs, pbs_here')
            push!(I_vals, I)
            if verbose
                println("Phase boundaries are at ", pbs_here)
            end
        catch e
            println("Solve failed at I=", I)
        end
    end
    return vcat(pbs[:,1], reverse(pbs[:,2])), vcat(I_vals, reverse(I_vals))
end
