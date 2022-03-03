using ElectrochemicalKinetics
using Plots
using NLsolve
# using QuadGK
using FastGaussQuadrature

# define some constants and default parameters
kB = 8.617e-5
T = 298
muoA = 0.02
muoB = 0.03
Ω = 0.1
N = 1000
quadfun = gausslegendre


# our familiar thermodynamic functions
hs(x;Ω=Ω) = @. x*(1-x)*Ω # enthalpy of mixing
s(x) = @. -kB*(x*log(x+eps(Float64)) + (1-x)*log(1-x+eps(Float64))) # entropy per particle...added epsilons in the hopes that things will be less obnoxious at the edges
g_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) = @. hs(x;Ω=Ω) - T*s(x)+ muoA*(1-x) + muoB*x # Gibbs free energy per particle
μ_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) = @. (1-2*x)*Ω + kB*T*log(x/(1-x)) + muoB-muoA # chemical potential

"""
µ and g with kinetic constributions, can be modeled using any <:KineticModel object

These functions return single-argument functions (to easily use common-tangent function below while
still being able to swap out model parameters by calling "function-builders" with different arguments).
"""
function µ_kinetic(I, km::KineticModel; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    thermo_term(x) = μ_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    μ(x::Real) = thermo_term(x) .+ fit_overpotential((1-x)*km, I)
    μ(x::AbstractVector) = thermo_term(x) .+ fit_overpotential((1 .- x).*Ref(km), I)
    return μ
end
function g_kinetic(I, km::KineticModel; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    thermo_term(x) = g_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    function kinetic_term(x, w, n)
        # f(x) = ElectrochemicalKinetics.fit_overpotential_t( (1 .- x) .* Ref(km), I, 0.1)
        f(x) = ElectrochemicalKinetics.fit_overpotential( (1 .- x) .* Ref(km), I, true)
        map((w, n) -> sum(w .* f(n)), eachcol(w), eachcol(n))
    end
    # g(x::AbstractVector) = thermo_term(x) .+ kinetic_term(x)
    g(x, w, n) = thermo_term(x) .+ kinetic_term(x, w, n)
    return g
end

# zeros of this function correspond to pairs of x's satisfying the common tangent condition for a given µ function
# case where we just feed in two points (x should be of length two)
function common_tangent(x::Vector, I, km::KineticModel; Ω=Ω, muoA=muoA, muoB=muoB, T=T, nodes, weights)
    g = g_kinetic(I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    µ = µ_kinetic(I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    [(g(x[2], weights, nodes) - g(x[1], weights, nodes))/(x[2] - x[1]) .- μ(x[1]), μ(x[2])-μ(x[1])]
end

# case where we want to check many points at once (shape of x should be N x 2)
function common_tangent(x::Array, I, km::KineticModel; nodes, weights, kwargs...)
    g = g_kinetic(I, km; kwargs...)
    µ = µ_kinetic(I, km; kwargs...)
    # w1, n1 = ElectrochemicalKinetics.scale_integration_nodes(zeros.(x[:, 1]), x[:, 1])
    # w2, n2 = ElectrochemicalKinetics.scale_integration_nodes(zero.(x[:, 2]), x[:, 2])
    w1, n1 = weights[:, 1], nodes[:, 1]
    w2, n2 = weights[:, 2], nodes[:, 2]
    Δg = g(x[:,2], w2, n2) - g(x[:,1], w1, n1)
    Δx = x[:,2] - x[:,1]
    μ1 = μ(x[:,1])
    Δμ = μ(x[:,2]) - μ1
    hcat(Δg ./ Δx .- μ1, Δμ)
end

function find_phase_boundaries(I, km::KineticModel; quadfun=quadfun, N=N, kwargs...)
    unscaled_x, unscaled_w = quadfun(N)
    
    function myct!(storage, x)
        nodes, weights = ElectrochemicalKinetics.scale_integration_nodes(unscaled_x,    # nodes
                                                                         unscaled_w,    # weights
                                                                         zero.(x),      # lb
                                                                         x)             # ub
        res = common_tangent(x, I, km; nodes = nodes, weights = weights, kwargs...)
        storage[1] = res[1]
        storage[2] = res[2]
    end
    # TODO: grad! here from AD
    x1 = nlsolve(myct!, [0.05 0.95])
    x1.zero
end

function phase_diagram(km::KineticModel; kwargs...)
    # TODO: write this, lol
end

# bv = ButlerVolmer(300)
# I_vals = 10 .^ (1.1:0.025:3.1)


# this line takes a few seconds with B-V but aaages with anything else...
# pbs = [find_phase_boundaries(I, bv, T=330) for I in I_vals]

# plot(vcat(pbs...), hcat(I_vals, I_vals), label="phase boundary")
# xlabel!("x")
# ylabel!("I")
