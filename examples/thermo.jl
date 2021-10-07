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

# enthalpy of mixing
@. hs(x;Ω=Ω) = x*(1-x)*Ω

# entropy per particle
@. s(x) = -kB*(x*log(x) + (1-x)*log(1-x))

# Gibbs free energy per particle
@. g_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) =
  hs(x;Ω=Ω) - T*s(x)+ muoA*(1-x) + muoB*x

# chemical potential
@. μ_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) =
  (1 - 2*x)*Ω + kB*T*log(x/(1-x)) + (muoB-muoA)

"""
µ and g with kinetic constributions, can be modeled using any <:KineticModel object

These functions return single-argument functions (to easily use common-tangent function below while
still being able to swap out model parameters by calling "function-builders" with different arguments).
"""
function µ_kinetic(I::Real, km::KineticModel;
                   Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    function (x)
      µ_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) .+
      first.(fit_overpotential.((1 .- x) .* Ref(km), I))
    end
end

function µ_kinetic(I::AbstractVector, km::KineticModel;
                   Ω=Ω, muoA=muoA, muoB=muoB, T=T)
  function (x)
    µ_thermo(x, Ω=Ω, muoA=muoA, muoB=muoB, T=T) .+
    first.(fit_overpotential.((1 .- x) .* Ref(km), Ref(I)))
  end
end

# g_kinetic(I, km::KineticModel; Ω=Ω, muoA=muoA, muoB=muoB, T=T) = 
#     x -> g_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) .+ [q[1][1][1] for q in quadgk.(y->fit_overpotential.((1 .-y).* repeat([km], length(y)), I), 0, x)]

# fully scalar version
function g_kinetic(I::Real, km::KineticModel;
                   Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    function (x)
      g_thermo(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T) +
      quadgk(y -> first(fit_overpotential((1 - y) * km, I)), 0, x)[1]
    end
end

function g_kinetic(I::AbstractVector, km::KineticModel;
                   Ω=Ω, muoA=muoA, muoB=muoB, T=T)
  g_kinetic.(I, Ref(km); Ω=Ω, muoA=muoA, muoB=muoB, T=T)
end

# for vector inputs x (for a vector of I's we would get a vector of functions that would then presumably each need to be broadcasted)
function g_kinetic_vecx(I::Real, km::KineticModel; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    function g_func(x::AbstractVector)
        @assert issorted(x) "Input x values must be sorted!"
        thermo_terms = g_thermo.(x; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
        lower_bounds = vcat([0], x[1:end-1])
        upper_bounds = x
        integral_slices = quadgk.(y -> fit_overpotential((1-y)*km, I), lower_bounds, upper_bounds)
        integral_slices_f = first.(integral_slices)
        kinetic_terms = cumsum([t[1] for t in integral_slices_f])
        thermo_terms .+ kinetic_terms
    end
end

function g_kinetic_vecxI(I::AbstractVector, km::KineticModel;
                         Ω=Ω, muoA=muoA, muoB=muoB, T=T)

  # @show "in g_kinetic_vecxI"
  function my_vec(x)
    # A = fill(km.A, length(I))
    # α = fill(km.α, length(I))
    thermo = g_thermo.(x; Ω, muoA, muoB, T)
    lower_bounds = vcat([0], x[1:end-1])
    upper_bounds = x
    @show length(I)
    integral_slices = quadgk.(y -> fit_overpotential((1-y)*km, I), lower_bounds, upper_bounds)
    @show typeof(integral_slices)
    integral_slices_f = first.(integral_slices)
    # @show integral_slices_f[1]
    # kinetic_terms = cumsum([t[1] for t in integral_slices_f])
    # @show kinetic_terms
    thermo .+ integral_slices_f[1] # kinetic_terms

  end
  # g_kinetic_vecx.(I, Ref(km); Ω=Ω, muoA=muoA, muoB=muoB, T=T)
  # my_vec(I) # , Ω=Ω, muoA=muoA, muoB=muoB, T=T)
end

# zeros of this function correspond to pairs of x's satisfying the common tangent condition for a given µ function
function common_tangent(x, I, km::KineticModel; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    g = g_kinetic_vecxI(I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    µ = µ_kinetic(I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    x1 = @view x[:, 1]
    x2 = @view x[:, 2]
    ((g(x2).-g(x1))./(x2.-x1) .- μ(x1))[1], (μ(x2).-μ(x1))[1]
end

# x should be N x 2, I of length N
function common_tangent_vecI(x, I::AbstractVector, km::KineticModel;
                             Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    # @show "in common_tangent_vecI"
    # gs = g_kinetic(I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    gs = g_kinetic_vecxI(I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    # @show "got gs"
    µs = µ_kinetic(I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    inds = 1:length(I)
    # r = x[:, 1]
    # @show gs(r)
    r1 = @views gs(x[:, 1])
    r2 = @views gs(x[:, 2])
    dg = r2 - r1 # map((g,r) -> g(r[2]) - g(r[1]), gs, eachrow(x))
    dx = @views x[:, 2] .- x[:, 1]
    μ1 = @views μs(x[:, 1])
    μ2 = @views μs(x[:, 2])
    # @show μ1
    # µ1 = map((m,r) -> m(r[1])[1], μs, eachrow(x))
    # µ2 = map((m,r) -> m(r[2])[1], μs, eachrow(x))
    return hcat(dg ./ dx .- µ1, µ2 .- µ1) # result N x 2
end

function find_phase_boundaries(I, km::KineticModel; Ω=Ω, muoA=muoA, muoB=muoB, T=T)
    # @show "in find_phase"
    function myct!(storage, x)
        # @show size(x)
        res = common_tangent_vecI(x, I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T) # N x 2
        storage .= res
    end
    # TODO: grad! here from AD
    function grad!(S, V)
      # @show size(S)
      # @show size(V)
      j = Zygote.jacobian(V -> common_tangent_vecI(V, I, km; Ω=Ω, muoA=muoA, muoB=muoB, T=T), V)
      # @show size(j[1])
      S .= j[1] # V
    end
    init = zeros(length(I), 2)
    init[:, 1] .= 0.05
    init[:, 2] .= 0.95
    x1 = nlsolve(myct!, grad!, init, show_trace = true)
    x1.zero
end

bv = ButlerVolmer(300)
# I_vals = 10 .^ (1.1:0.025:3.1)
I_vals = 10 .^ (1.1:0.1:3.1)

init = zeros(length(I_vals), 2)
init[:, 1] .= 0.05
init[:, 2] .= 0.95

# this line takes a few seconds with B-V but aaages with anything else...
# pbs = [find_phase_boundaries(I, bv, T=330) for I in I_vals]
# pbs = find_phase_boundaries(I_vals, bv, T=330)

#plot(vcat(pbs...), hcat(I_vals, I_vals), label="phase boundary")
# plot(pbs, hcat(I_vals, I_vals), label="phase boundary")
# xlabel!("x")
# ylabel!("I")
