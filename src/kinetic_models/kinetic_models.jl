using .DOS
using SpecialFunctions
using Statistics
using Interpolations
using DelimitedFiles

"""
    fermi_dirac(E, kT=0.026)

Compute the value of the Fermi-Dirac distribution at energy `E` (relative to the Fermi energy) and thermal energy `kT`.
"""
fermi_dirac(E; kT = 0.026) = inv.(1 .+ exp.(E ./ kT))

"""
    KineticModel

Abstract base class for kinetic models. All subtypes need to dispatch the `compute_k` function.
"""
abstract type KineticModel end

# dispatch for net rates, returns absolute value
(km::KineticModel)(V_app; kT = 0.026) =
    abs.(km(V_app, Val(true); kT = kT) - km(V_app, Val(false); kT = kT))

# return a new one with a scaled prefactor
import Base.*
function *(c::Real, km::KineticModel)
    new_A = c*km.A
    other_params = getfield.([km], propertynames(km))[2:end]
    typeof(km)(new_A, other_params...)
end

# generic pretty printing
function Base.show(io::IO, m::KineticModel)
    s = repr(typeof(m)) * "("
    for field in propertynames(m)
        s =
            s *
            string(field) *
            "=" *
            string(round(getproperty(m, field), sigdigits = 3)) *
            ", "
    end
    s = s[1:end-2] * ")"
    print(io, s)
end

"""
    IntegralModel

Abstract base class for "Marcus-like" kinetic models that require computation of an energy integral. All subtypes need to dispatch the `compute_k` function directly, or dispatch the `integrand` function and use the default `compute_k` dispatch.
"""
abstract type IntegralModel <: KineticModel end # "Marcus-like"

# check to catch missed dispatches for new types
# integrand(km::IntegralModel, V_dl, ox::Bool; kwargs...) =
#     error("An integral-based kinetic model must dispatch the `integrand` function!")

# TODO: check that this passes through both kT and V_q appropriately
# dispatch for net rates
integrand(km::IntegralModel, V_dl; kwargs...) =
    E -> abs.(
        integrand(km, V_dl, true; kwargs...)(E) - integrand(km, V_dl, false; kwargs...)(E),
    )

    """
    compute_k(V_app, model::KineticModel, ox::Bool; kwargs...)
    compute_k(V_app, model::KineticModel; kwargs...)
    compute_k(E_min, E_max, V_app, model::IntegralModel, ox::Bool; kwargs...)
    compute_k(E_min, E_max, V_app, model::IntegralModel; kwargs...)
    compute_k(E_min, E_max, V_app, model::MarcusHushChidseyDOS, calc_cq::Bool=false; C_dl = 10.0, Vq_min = -0.5, Vq_max = 0.5, kwargs...)

Compute the rate constant k predicted by a given kinetic model at a applied voltage `V_app`. If a flag for reaction direction `ox` is supplied, `true` gives the oxidative and `false` the reductive direction, while omitting this flag will yield net reaction rate (absolute value thereof).

If the model is an `IntegralModel`, integration bounds `E_min` and `E_max` must be supplied. Integration is done via GK quadrature.

If calc_cq flag is passed, optionally compute voltage shifts due to quantum capacitance.
"""
compute_k(V_app, model::KineticModel, args...; kT = 0.026) = model(V_app, args...; kT = kT)

compute_k(
    V_app,
    model::IntegralModel,
    args...; # would just be the ox flag
    kT = 0.026,
    E_min = -100 * kT,
    E_max = 100 * kT,
    kwargs...,
) = begin
  n, w = scale(E_min, E_max)
  f = integrand(model, V_app, args...; kT = kT)
  sum(w .* f.(n))
end

include("ButlerVolmer.jl")
include("Marcus.jl")
include("AsymptoticMarcusHushChidsey.jl")
include("MarcusHushChidsey.jl")
include("MarcusHushChidseyDOS.jl")


