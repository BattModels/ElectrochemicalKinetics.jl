using .DOS

"""
    fermi_dirac(E, kT=0.026)

Compute the value of the Fermi-Dirac distribution at energy `E` (relative to the Fermi energy) and thermal energy `kT`.
"""
fermi_dirac(E; kT = 0.026) = 1 / (1 + exp(E / kT))

abstract type KineticModel end

# check to catch missed dispatches for new types
integrand(km::KineticModel, V_dl::Real, ox::Bool; kT::Real = 0.026) =
    error("A kinetic model must dispatch the `integrand` function!")

# TODO: check that this passes through both kT and V_q appropriately
# dispatch for net rates
integrand(km::KineticModel, V_dl::Real; args...) =
    E -> abs(integrand(km, V_dl, true; args...)(E) - integrand(km, V_dl, false; args...)(E))

struct Marcus <: KineticModel
    λ::Float64
end

function integrand(m::Marcus, V_dl::Real, ox::Bool; kT::Real = 0.026)
    arg =
        E ->
            ox ? -(m.λ - V_dl + E)^2 / (4 * m.λ * kT) : -(m.λ + V_dl - E)^2 / (4 * m.λ * kT)
    E -> exp(arg(E))
end

struct MarcusHushChidsey <: KineticModel
    λ::Float64
    average_dos::Float64
end

MarcusHushChidsey(λ, dd::DOSData) = MarcusHushChidsey(λ, dd.average_value)
MarcusHushChidsey(λ, dos_file::String; args...) =
    MarcusHushChidsey(λ, DOSData(dos_file; args...))

function integrand(mhc::MarcusHushChidsey, V_dl::Real, ox::Bool; kT::Real = 0.026)
    marcus = integrand(Marcus(mhc.λ), V_dl, ox; kT = kT)
    fd(E) = ox ? 1 - fermi_dirac(E; kT = kT) : fermi_dirac(E; kT = kT)
    return E -> mhc.average_dos * marcus(E) * fd(E)
end

struct MarcusHushChidseyDOS <: KineticModel
    λ::Float64
    dos::DOSData
end

MarcusHushChidseyDOS(λ, dos_file::String; args...) =
    MarcusHushChidseyDOS(λ, DOSData(dos_file; args...))

integrand(mhcd::MarcusHushChidseyDOS, V_dl::Real, ox::Bool; kT::Real = 0.026, V_q = 0.0) =
    E ->
        mhcd.dos.interp_func(E + V_q) *
        integrand(MarcusHushChidsey(mhcd.λ, 1.0), V_dl, ox; kT = kT)(E)
