using .DOS
using SpecialFunctions
using Statistics
using Interpolations
using DelimitedFiles
using TimerOutputs

const to = TimerOutput()

include("../utils/misc.jl")

"""
    fermi_dirac(E, T=298)

Compute the value of the Fermi-Dirac distribution at energy `E` (relative to the Fermi energy) and temperature `T`.
"""
fermi_dirac(E; T = 298) = inv.(1 .+ exp.(E ./ (kB * T)))

"""
    KineticModel

Abstract base class for kinetic models.
"""
abstract type KineticModel{T} end

# return a new one with a scaled prefactor
import Base.*
function *(c, km::KineticModel)
    new_A = c*km.A
    other_params = getfield.([km], propertynames(km))[2:end]
    eval(nameof(typeof(km)))(new_A, other_params...)
end

import Base.length
length(km::KineticModel) = length(km.A)

# generic pretty printing
function Base.show(io::IO, m::KineticModel)
    s = repr(typeof(m)) * "("
    for field in propertynames(m)
        s =
            s *
            string(field) *
            "=" *
            string(round.(getproperty(m, field), sigdigits = 3)) *
            ", "
    end
    s = s[1:end-2] * ")"
    print(io, s)
end

"""
    NonIntegralModel

Abstract base class for kinetic models whose rates can be computed directly from an input voltage without requiring an energy integral. All subtypes must dispatch the `rate_constant` function.
"""
abstract type NonIntegralModel{T} <: KineticModel{T} end

"""
    IntegralModel

Abstract base class for "Marcus-like" kinetic models that require computation of an energy integral. All subtypes need to dispatch the `rate_constant` function directly, or dispatch the `integrand` function and use the default `rate_constant` dispatch.
"""
abstract type IntegralModel{T} <: KineticModel{T} end

# check to catch missed dispatches for new types
# integrand(km::IntegralModel, V_dl, ox::Bool; kwargs...) =
#     error("An integral-based kinetic model must dispatch the `integrand` function!")

# TODO: check that this passes through V_q appropriately
# dispatch for net rates
function integrand(km::IntegralModel, V; a_r=1.0, a_o=1.0, T=298, kwargs...)
    @timeit to "net rate integrand" begin
        eq_V = kB * T .* log.(a_r ./ a_o)
        E -> a_r .* integrand(km, V .- eq_V, Val(true); kwargs...)(E) .- a_o .* integrand(km, V .- eq_V, Val(false); kwargs...)(E)
    end
end

integrand(km::IntegralModel, V, ox::Bool; kwargs...) = integrand(km, V, Val(ox); kwargs...)

"""
    rate_constant(V_app, model::KineticModel, ox::Bool; kwargs...)
    rate_constant(V_app, model::KineticModel; a_r=1.0, a_o=1.0, kwargs...)
    rate_constant(E_min, E_max, V_app, model::MarcusHushChidseyDOS, calc_cq::Bool=false; C_dl = 10.0, Vq_min = -0.5, Vq_max = 0.5, kwargs...)

Compute the rate constant k predicted by a given kinetic model at a applied voltage `V_app`. If a flag for reaction direction `ox` is supplied, `true` gives the oxidative and `false` the reductive direction, while omitting this flag will yield net reaction rate. In the net rate case, activities of reduced and oxidized species (`a_r` and `a_o`, respectively) may also be supplied (default values are unity).

If the model is an `IntegralModel`, integration bounds `E_min` and `E_max` may be supplied as kwargs. Integration is done via GK quadrature.

If calc_cq flag is passed, optionally compute voltage shifts due to quantum capacitance (only applicable to `MarcusHushChidseyDOS` models).
"""
rate_constant(V_app, model::NonIntegralModel, ox::Bool; kwargs...) = rate_constant(V_app, model, Val(ox); kwargs...)

# dispatch for net rates
function rate_constant(V_app, model::NonIntegralModel; a_r=1.0, a_o=1.0, T=298, kwargs...)
    @timeit to "net rate nonintegral" begin
        eq_V = kB * T .* log.(a_r./a_o)
        a_r .* rate_constant(V_app .- eq_V, model, Val(true); T=T, kwargs...) .- a_o .* rate_constant(V_app .- eq_V, model, Val(false); T=T, kwargs...)
    end
end

# TODO: add tests that both args and kwargs are correctly captured here (also for the Val thing)
# "callable" syntax
(m::KineticModel)(V_app, args...; kwargs...) = rate_constant(V_app, m, args...; kwargs...)

function rate_constant(
    V_app,
    model::IntegralModel,
    args...; # would just be the ox flag, if present
    T = 298,
    E_min = -100 * kB * T,
    E_max = 100 * kB * T,
    kwargs... # e.g. a_o, a_r
)
    n, w = scale(E_min, E_max)
    f = integrand(model, V_app, args...; T = T, kwargs...)
    sum(w .* f.(n))
end

include("ButlerVolmer.jl")
include("Marcus.jl")
include("AsymptoticMarcusHushChidsey.jl")
include("MarcusHushChidsey.jl")
include("MarcusHushChidseyDOS.jl")
