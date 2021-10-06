using ElectrochemicalKinetics
using Plots
using NLsolve
using QuadGK

# define some constants and default parameters
kB = 8.617e-5
T = 298
muoA = 0.02
muoB = 0.03
Ω = 0.1

# our familiar thermodynamic functions
@. hs(x;Ω=Ω) = x*(1-x)*Ω # enthalpy of mixing
@. s(x) = -kB*(x*log(x) + (1-x)*log(1-x)) # entropy per particle
@. g_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) = hs(x;Ω=Ω) - T*s(x)+ muoA*(1-x) + muoB*x # Gibbs free energy per particle
@. μ_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) = (1-2*x)*Ω + kB*T*log(x/(1-x)) + muoB-muoA # chemical potential

"""
µ and g with kinetic constributions, can be modeled using any <:KineticModel object

These functions return single-argument functions (to easily use common-tangent function below while
still being able to swap out model parameters by calling "function-builders" with different arguments).
"""
µ_kinetic(I, km::KineticModel; Ω=Ω, muoA=muoA, muoB=muoB, T=T) = 
    x -> µ_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) .+ [t[1] for t in fit_overpotential.((1 .-x).* repeat([km], length(x)), I)]
# g_kinetic(I, km::KineticModel; Ω=Ω, muoA=muoA, muoB=muoB, T=T) = 
#     x -> g_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) .+ [q[1][1][1] for q in quadgk.(y->fit_overpotential.((1 .-y).* repeat([km], length(y)), I), 0, x)]

# fully scalar version
g_kinetic(I, km::KineticModel; Ω=Ω, muoA=muoA, muoB=muoB, T=T) = 
    x -> g_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) + quadgk(y->fit_overpotential((1 - y) * km, I), 0, x)[1][1]

# for vector inputs x (for a vector of I's we would get a vector of functions that would then presumably each need to be broadcasted)
function g_kinetic_vec(I, km::KineticModel; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    function g_func(x::Vector{<:Real})
        @assert issorted(x) "Input x values must be sorted!"
        thermo_terms = g_thermo.(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
        lower_bounds = vcat([0], x[1:end-1])
        upper_bounds = x
        integral_slices = quadgk.(y -> fit_overpotential((1-y)*km, I), lower_bounds, upper_bounds)
        kinetic_terms = cumsum([t[1][1] for t in integral_slices])
        thermo_terms .+ kinetic_terms
    end
    return g_func
end

# zeros of this function correspond to pairs of x's satisfying the common tangent condition for a given µ function
function common_tangent(x, I, km::KineticModel; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    g = g_kinetic(I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    µ = µ_kinetic(I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    [((g(x[:,2]).-g(x[:,1]))./(x[:,2].-x[:,1]) .- μ(x[:,1]))[1], (μ(x[:,2]).-μ(x[:,1]))[1]]
end

function find_phase_boundaries(I, km::KineticModel; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    function myct!(storage, x)
        res = common_tangent(x, I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
        storage[1] = res[1]
        storage[2] = res[2]
    end
    # TODO: grad! here from AD
    x1 = nlsolve(myct!, [0.05 0.95])
    x1.zero
end


bv = ButlerVolmer(300)
I_vals = 10 .^ (1.1:0.025:3.1)

# this line takes a few seconds with B-V but aaages with anything else...
pbs = [find_phase_boundaries(I, bv, T=330) for I in I_vals]

plot(vcat(pbs...), hcat(I_vals, I_vals), label="phase boundary")
xlabel!("x")
ylabel!("I")