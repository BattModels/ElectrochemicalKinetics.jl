using ElectrochemicalKinetics
using Plots
using NLsolve
using FastGaussQuadrature

# define some constants and default parameters
const room_T = 298
const muoA_default = 0.02
const muoB_default = 0.03
const Ω_default = 0.1
N = 1000
quadfun = gausslegendre

# our familiar thermodynamic functions
h(x;Ω=Ω_default) = @. x*(1-x)*Ω # enthalpy of mixing
s(x) = @. -kB*(x*log(x+eps(Float64)) + (1-x)*log(1-x+eps(Float64))) # entropy per particle...added epsilons in the hopes that things will be less obnoxious at the edges
g_thermo(x; Ω=Ω_default, muoA=muoA_default, muoB=muoB_default, T=room_T) = @. h(x;Ω=Ω) - T*s(x)+ muoA*(1-x) + muoB*x # Gibbs free energy per particle
μ_thermo(x; Ω=Ω_default, muoA=muoA_default, muoB=muoB_default, T=room_T) = @. (1-2*x)*Ω + kB*T*log(x/(1-x)) + muoB-muoA # chemical potential

"""
µ and g with kinetic constributions, can be modeled using any <:KineticModel object
"""
function µ_kinetic(x, I, km::KineticModel; intercalate=true, kwargs...)
    prefactor(x) = intercalate ? (1 .- x) : x
    return μ_thermo(x; kwargs...) .+ overpotential(I, km, prefactor(x))
end

# TODO: test cases of vector/scalar x/I
function g_kinetic(x, I, km::KineticModel; intercalate=true, stepsize=1e-3, kwargs...)
    prefactor(x) = intercalate ? (1 .- x) : x
    # g(x) = thermo_term(x) .+ kinetic_term(vec(x))
    # g(x::Real) = thermo_term(x) + kinetic_term(x)[1]
    f(x) = ElectrochemicalKinetics.overpotential(I, km, prefactor(x); kwargs...)
    # n, w = ElectrochemicalKinetics.scale(zero.(x), x)
    n = 0:stepsize:maximum(x)
    f_vals = f(n)

    # trapezoidal rule
    # f_vals = 0.5 .* (f_vals[1:end-1] .+ f_vals[2:end])
    # n = n[2:end]

    # Simpson's 3/8
    f_vals = (1/8) .* (f_vals[1:end-3] .+ 3*f_vals[2:end-2] .+ 3*f_vals[3:end-1] .+ f_vals[4:end])
    n = n[4:end]

    w = stepsize
    stops = map(xv->searchsorted(n, xv).stop, x)
    kinetic_part = w .* map(x->sum(f_vals[1:x]), stops)
    return g_thermo(x; kwargs...) + kinetic_part
    # return g_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) .+ map((w, n) -> sum(w .* f(n)), eachcol(w), eachcol(n))
    # return map((w, n) -> sum(w .* f(n)), eachcol(w), eachcol(n)) # just kinetic part for now
end

#TODO: finish this, need term for thermo part in derivative wrt x and also derivative wrt I, currently set to nothing which is very much incorrect
# Zygote.@adjoint function g_kinetic(x, I, km::KineticModel; kwargs...)
#     y, pb = Zygote.pullback(x -> g_kinetic(x, I, km::KineticModel; kwargs...), I) # I want "default" behavior wrt I
#     y, function(D) # I don't understand what `D` represents above so I'm just copying it
#         custom_grad = value_at_bound(x, I, km; kwargs...) # this is my Leibniz rule thing
#         (custom_grad, pb(D), nothing) # x, I, km args give adjoints of custom, default, nothing
#     end
# end

# zeros of this function correspond to pairs of x's satisfying the common tangent condition for a given µ function
# case where we just feed in two points (x should be of length two)
# probably don't actually need this separate version, keeping it here for sanity check for now
function common_tangent(x::Vector, I, km::KineticModel; stepsize=5e-3, kwargs...)
    g = g_kinetic(x, I, km; stepsize=stepsize, kwargs...)
    µ = µ_kinetic(x, I, km; kwargs...)
    [(g[2] - g[1])/(x[2] - x[1]) - μ[1], μ[2]-μ[1]]
end

# case where we want to check many points at once (shape of x should be N x 2)
# TODO: check that this works properly
function common_tangent(x::Array, I, km::KineticModel; kwargs...)
    g(x) = g_kinetic(x, I, km; kwargs...)
    µ(x) = µ_kinetic(x, I, km; kwargs...)
    Δg = g(x[:,2]) .- g(x[:,1])
    Δx = x[:,2] .- x[:,1]
    μ1 = μ(x[:,1])
    Δμ = μ(x[:,2]) .- μ1
    hcat(Δg ./ Δx .- μ1, Δμ)
end

# TODO: get vector I case working
# also even scalar is super slow right now
function find_phase_boundaries(I::Real, km::KineticModel; stepsize=5e-3, verbose=false, kwargs...)
    
    function myct!(storage, x)
        res = common_tangent(x, I, km; stepsize=stepsize, kwargs...)
        storage[1] = res[1]
        storage[2] = res[2]
    end
    # TODO: grad! here from AD
    x1 = nlsolve(myct!, [0.05, 0.95], show_trace=verbose)
    x1.zero
end

function phase_diagram(km::KineticModel; kwargs...)
    # TODO: write this, lol
end

# bv = ButlerVolmer(300, 0.5)
# I_vals = 10 .^ (1.1:0.1:3.1)


# this line takes a few seconds with B-V but aaages with anything else...
# pbs = [find_phase_boundaries(I, bv, T=330) for I in I_vals]

# plot(vcat(pbs...), hcat(I_vals, I_vals), label="phase boundary")
# xlabel!("x")
# ylabel!("I")
