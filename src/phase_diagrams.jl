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

"""
µ and g with kinetic constributions, can be modeled using any <:KineticModel object

These functions return single-argument functions (to easily use common-tangent function below while
still being able to swap out model parameters by calling "function-builders" with different arguments).
"""
function µ_kinetic(I, km::KineticModel; Ω=Ω_default, muoA=muoA_default, muoB=muoB_default, T=room_T)
    thermo_term(x) = μ_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    μ(x::Real) = thermo_term(x) .+ overpotential(I, (1-x)*km)
    μ(x::AbstractVector) = thermo_term(x) .+ overpotential(I, (1 .- x).*Ref(km))
    return μ
end

function g_kinetic(I, km::KineticModel; Ω=Ω_default, muoA=muoA_default, muoB=muoB_default, T=room_T)
    thermo_term(x) = g_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    #TODO: gradient of this term is just value of overpotential(x)
    function kinetic_term(x)
        f(x) = ElectrochemicalKinetics.overpotential(I, (1 .- x) .* Ref(km))
        n, w = ElectrochemicalKinetics.scale_coarse(zero.(x), x)
        map((w, n) -> sum(w .* f(n)), eachcol(w), eachcol(n))
    end
    g(x) = thermo_term(x) .+ kinetic_term(vec(x))
    g(x::Real) = thermo_term(x) + kinetic_term(x)[1]
    return g
end

# zeros of this function correspond to pairs of x's satisfying the common tangent condition for a given µ function
# case where we just feed in two points (x should be of length two)
function common_tangent(x::Vector, I, km::KineticModel; Ω=Ω_default, muoA=muoA_default, muoB=muoB_default, T=room_T)
    g = g_kinetic(I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    µ = µ_kinetic(I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    [(g(x[2]) - g(x[1]))/(x[2] - x[1]) .- μ(x[1]), μ(x[2])-μ(x[1])]
end

# case where we want to check many points at once (shape of x should be N x 2)
# TODO: test this...may not really be able to save much time
function common_tangent(x::Array, I, km::KineticModel; kwargs...)
    g = g_kinetic(I, km; kwargs...)
    µ = µ_kinetic(I, km; kwargs...)
    Δg = g(x[:,2]) - g(x[:,1])
    Δx = x[:,2] - x[:,1]
    μ1 = μ(x[:,1])
    Δμ = μ(x[:,2]) - μ1
    hcat(Δg ./ Δx .- μ1, Δμ)
end

function find_phase_boundaries(I, km::KineticModel; guess=[0.05, 0.95], verbose=false, kwargs...)
    
    function myct!(storage, x)
        res = common_tangent(x, I, km; kwargs...)
        storage[1] = res[1]
        storage[2] = res[2]
    end
    # TODO: grad! here from AD
    x1 = nlsolve(myct!, guess, show_trace=verbose, ftol=1e-6, xtol=1e-6)
    x1.zero
end

function phase_diagram(km::KineticModel; kwargs...)
    # TODO: write this, lol
end
