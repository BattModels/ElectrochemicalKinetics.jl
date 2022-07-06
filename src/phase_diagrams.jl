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
    f(x) = overpotential(I, km, prefactor(x); kwargs...)
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
    return g_thermo(x; kwargs...) + w .* map(x->sum(f_vals[1:x]), stops)
    # return g_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) .+ map((w, n) -> sum(w .* f(n)), eachcol(w), eachcol(n))
    # return map((w, n) -> sum(w .* f(n)), eachcol(w), eachcol(n)) # just kinetic part for now
end

#TODO: finish/test this
Zygote.@adjoint function g_kinetic(x, I, km::KineticModel; intercalate=true, stepsize=1e-3, kwargs...)
    y, pb = Zygote.pullback(j -> g_kinetic(x, j, km; intercalate=intercalate, stepsize=stepsize, kwargs...), I) # "default" behavior wrt I
    y, function(Δ)
        prefactor(x) = intercalate ? (1 .- x) : x
        custom_grad = overpotential(I, km, prefactor(x); kwargs...) # Leibniz rule
        (Δ .* custom_grad, Δ .* pb(Δ), nothing)
    end
end

# zeros of this function correspond to pairs of x's satisfying the common tangent condition for a given µ function
# case where we just feed in two points (x should be of length two)
# probably don't actually need this separate version, keeping it here for sanity check for now
function common_tangent(x::Vector, I, km::KineticModel; stepsize=5e-3, kwargs...)
    @timeit to "common_tangent (vector)" begin
    @timeit to "Δg" g = g_kinetic(x, I, km; stepsize=stepsize, kwargs...)
    @timeit to "Δµ" µ = µ_kinetic(x, I, km; kwargs...)
    [(g[2] - g[1])/(x[2] - x[1]) - μ[1], μ[2]-μ[1]]
    end
end

# case where we want to check many points at once (shape of x should be N x 2)
# TODO: check that this works properly
function common_tangent(x::Array, I, km::KineticModel; kwargs...)
    @timeit to "common_tangent (array)" begin
    g(x) = g_kinetic(x, I, km; kwargs...)
    µ(x) = µ_kinetic(x, I, km; kwargs...)
    @timeit to "Δg" Δg = g(x[:,2]) .- g(x[:,1])
    Δx = x[:,2] .- x[:,1]
    @timeit to "Δµ" begin
    μ1 = μ(x[:,1])
    Δμ = μ(x[:,2]) .- μ1
    end
    hcat(Δg ./ Δx .- μ1, Δμ)
    end
end

# TODO: get vector I case working
# also even scalar is super slow right now
function find_phase_boundaries(I::Real, km::KineticModel; stepsize=5e-3, verbose=false, guess=[0.05, 0.95], kwargs...)
    @timeit to "find_phase_boundaries" begin
    function myct!(storage, x)
        res = common_tangent(x, I, km; stepsize=stepsize, kwargs...)
        storage[1] = res[1]
        storage[2] = res[2]
    end
    # TODO: grad! here from AD
    x1 = nlsolve(myct!, guess, show_trace=verbose)
    x1.zero
    end
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
