using .DOS
using SpecialFunctions

"""
    fermi_dirac(E, kT=0.026)

Compute the value of the Fermi-Dirac distribution at energy `E` (relative to the Fermi energy) and thermal energy `kT`.
"""
fermi_dirac(E; kT = 0.026) = 1 / (1 + exp(E / kT))

abstract type KineticModel end

# dispatch for net rates
(km::KineticModel)(V_app; kT = 0.026) =
    abs(km(V_app, true; kT = kT) - km(V_app, false; kT = kT))

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
    ButlerVolmer(A, α)
    ButlerVolmer(A)

Computes Butler-Volmer kinetics. 

If initialized with one argument, assumes symmetric electron transfer (α=0.5) and sets this to be the prefactor A. Note that this prefactor implicitly contains information about equilibrium activation energies, as well as geometric information.
"""
struct ButlerVolmer <: KineticModel
    A::Float64
    α::Float64
end

# default to symmetric response
ButlerVolmer(A) = ButlerVolmer(A, 0.5)

function (bv::ButlerVolmer)(V_app, ox::Bool; kT::Real = 0.026)
    local exp_arg
    if ox
        exp_arg = (bv.α * V_app) / kT
    else
        exp_arg = -((1 - bv.α) * V_app) / kT
    end
    bv.A * exp(exp_arg)
end

"""
    AsymptoticMarcusHushChidsey(A, λ)
    AsymptoticMarchusHushChidsey(λ)

Computes asymptotic solution to MHC model, as described in Zeng et al.: 10.1016/j.jelechem.2014.09.038d

If initialized with one argument, assumes this to be the reorganization energy λ and sets the prefactor to 1.0.
"""
struct AsymptoticMarcusHushChidsey <: KineticModel
    A::Float64
    λ::Float64
end

# default prefactor is 1
AsymptoticMarcusHushChidsey(λ) = AsymptoticMarcusHushChidsey(1.0, λ)

function (amhc::AsymptoticMarcusHushChidsey)(V_app, ox::Bool; kT::Real = 0.026)
    a = 1 + sqrt(amhc.λ)
    η = (2 * ox - 1) * V_app / kT
    λ_nondim = amhc.λ / kT
    arg = (λ_nondim - sqrt(a + η^2)) / (2 * sqrt(λ_nondim))
    pref = sqrt(π * λ_nondim) / (1 + exp(-η))
    return amhc.A * pref * erfc(arg)
end

abstract type IntegralModel <: KineticModel end # "Marcus-like"

# check to catch missed dispatches for new types
integrand(km::IntegralModel, V_dl::Real, ox::Bool; kwargs...) =
    error("An integral-based kinetic model must dispatch the `integrand` function!")

# TODO: check that this passes through both kT and V_q appropriately
# dispatch for net rates
integrand(km::IntegralModel, V_dl::Real; kwargs...) =
    E -> abs(
        integrand(km, V_dl, true; kwargs...)(E) - integrand(km, V_dl, false; kwargs...)(E),
    )

"""
    Marcus(A, λ)
    Marcus(λ)

Computes Marcus kinetics.

If initialized with one argument, assumes this to be the reorganization energy λ and sets the prefactor A to 1.0.
"""
struct Marcus <: KineticModel
    A::Float64
    λ::Float64
end

# default prefactor is 1
Marcus(λ) = Marcus(1.0, λ)

function (m::Marcus)(V_app, ox::Bool; kT::Real = 0.026)
    local exp_arg
    if ox
        exp_arg = -(m.λ + V_app)^2 / (4 * m.λ * kT)
    else
        exp_arg = -(m.λ - V_app)^2 / (4 * m.λ * kT)
    end
    m.A * exp(exp_arg)
end

"""
    MarcusHushChidsey(A, λ, average_dos)
    MarcusHushChidsey(λ, average_dos)

Computes Marcus-Hush-Chidsey kinetics: 10.1126/science.251.4996.919

Note that strictly speaking, `average_dos` and the prefactor `A` are redundant. They are both included primarily to facilitate comparisons with similarly parametrized `Marcus` models.
"""
struct MarcusHushChidsey <: IntegralModel
    A::Float64
    λ::Float64
    average_dos::Float64
end

# default prefactor is 1
MarcusHushChidsey(λ, avg_dos) = MarcusHushChidsey(1.0, λ, avg_dos)
# convert more detailed DOS information to just pull out average
MarcusHushChidsey(A, λ, dd::DOSData) = MarcusHushChidsey(A, λ, dd.average_value)
MarcusHushChidsey(A, λ, dos_file::String; kwargs...) =
    MarcusHushChidsey(A, λ, DOSData(dos_file; kwargs...))

# TODO: Check that both this and +DOS versions still match original paper
function integrand(mhc::MarcusHushChidsey, V_dl::Real, ox::Bool; kT::Real = 0.026)
    function marcus_term(E)
        local exp_arg
        if ox
            exp_arg = -(E-mhc.λ+V_dl)^2 / (4 * mhc.λ * kT)
        else
            exp_arg = -(E-mhc.λ-V_dl)^2 / (4 * mhc.λ * kT)
        end
        exp(exp_arg)
    end
    E -> mhc.A * mhc.average_dos * marcus_term(E) * fermi_dirac(E; kT = kT)
end

"""
    MarcusHushChidseyDOS(A, λ, dos)

Computes Marcus-Hush-Chidsey + DOS kinetics as described in Kurchin and Viswanathan: 10.1063/5.0023611 
"""
struct MarcusHushChidseyDOS <: IntegralModel
    A::Float64
    λ::Float64
    dos::DOSData
end

function Base.show(io::IO, mhcd::MarcusHushChidseyDOS)
    s = repr(typeof(mhcd)) * "("
    s *= "A=$(round(mhcd.A, sigdigits=3)), λ=$(round(mhcd.λ, sigdigits=3)))"
    print(io, s)
end

# default prefactor is 1
MarcusHushChidseyDOS(λ, dd::DOSData) = MarcusHushChidseyDOS(1.0, λ, dd)

MarcusHushChidseyDOS(A, λ, dos_file::String; kwargs...) =
    MarcusHushChidseyDOS(A, λ, DOSData(dos_file; kwargs...))

function integrand(
    mhcd::MarcusHushChidseyDOS,
    V_dl::Real,
    ox::Bool;
    kT::Real = 0.026,
    V_q = 0.0,
)
    function marcus_term(E)
        local exp_arg
        if ox
            exp_arg = -(mhcd.λ-V_dl+E)^2 / (4 * mhcd.λ * kT)
        else
            exp_arg = -(mhcd.λ+V_dl-E)^2 / (4 * mhcd.λ * kT)
        end
        exp(exp_arg)
    end
    fd(E) = ox ? 1 - fermi_dirac(E; kT = kT) : fermi_dirac(E; kT = kT)
    E -> mhcd.A * mhcd.dos.interp_func(E + V_q) * fd(E) * marcus_term(E)
end
