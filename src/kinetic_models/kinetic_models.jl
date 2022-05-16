using .DOS
using SpecialFunctions

"""
    fermi_dirac(E, kT=0.026)

Compute the value of the Fermi-Dirac distribution at energy `E` (relative to the Fermi energy) and thermal energy `kT`.
"""
fermi_dirac(E; kT = 0.026) = inv.(1 .+ exp.(E ./ kT))

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


include("ButlerVolmer.jl")
include("Marcus.jl")
include("AsymptoticMarcusHushChidsey.jl")
include("MarcusHushChidsey.jl")
include("MarcusHushChidseyDOS.jl")


