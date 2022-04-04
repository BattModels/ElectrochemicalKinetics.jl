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

# default to unit prefactor and symmetric response
ButlerVolmer() = ButlerVolmer(1.0, 0.5)
ButlerVolmer(A) = ButlerVolmer(A, 0.5)

function (bv::ButlerVolmer)(V_app, ::Val{true}; kT::Real = 0.026)
    exp_arg = (bv.α .* V_app) ./ kT
    bv.A .* exp.(exp_arg)
end

function (bv::ButlerVolmer)(V_app, ::Val{false}; kT::Real = 0.026)
    exp_arg = -((1 - bv.α) .* V_app) ./ kT
    bv.A .* exp.(exp_arg)
end

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
        exp_arg = -(m.λ .+ V_app).^2 ./ (4 * m.λ * kT)
    else
        exp_arg = -(m.λ .- V_app).^2 ./ (4 * m.λ * kT)
    end
    m.A .* exp.(exp_arg)
end

function (m::Marcus)(V_app, ::Val{true}; kT::Real = 0.026)
    exp_arg = -(m.λ .+ V_app).^2 ./ (4 * m.λ * kT)
    m.A .* exp.(exp_arg)
end

function (m::Marcus)(V_app, ::Val{false}; kT::Real = 0.026)
    exp_arg = -(m.λ .- V_app).^2 ./ (4 * m.λ * kT)
    m.A .* exp.(exp_arg)
end

"""
    AsymptoticMarcusHushChidsey(A, λ)
    AsymptoticMarcusHushChidsey(λ)

Computes asymptotic solution to MHC model, as described in Zeng et al.: 10.1016/j.jelechem.2014.09.038d, with a corrction prefactor of kT since there is an error in the nondimensionalization in that work.

If initialized with one argument, assumes this to be the reorganization energy λ and sets the prefactor to 1.0.
"""
struct AsymptoticMarcusHushChidsey <: KineticModel
    A::Float64
    λ::Float64
end

# default prefactor is 1
AsymptoticMarcusHushChidsey(λ) = AsymptoticMarcusHushChidsey(1.0, λ)

function (amhc::AsymptoticMarcusHushChidsey)(V_app, ox::Bool; kT::Real = 0.026)
    η = (2 * ox - 1) .* V_app ./ kT
    λ_nondim = amhc.λ / kT
    a = 1 + sqrt(λ_nondim)
    arg = (λ_nondim .- sqrt.(a .+ η.^2)) ./ (2 * sqrt(λ_nondim))
    pref = sqrt(π * λ_nondim) ./ (1 .+ exp.(-η))
    return kT * amhc.A .* pref .* erfc.(arg)
end


function (amhc::AsymptoticMarcusHushChidsey)(V_app; kT::Real = 0.026)
    η = V_app ./ kT
    λ_nondim = amhc.λ / kT
    a = 1 + sqrt(λ_nondim)
    arg = (λ_nondim .- sqrt.(a .+ η.^2)) ./ (2 * sqrt(λ_nondim))
    pref = sqrt(π * λ_nondim) .* tanh.(η/2)
    return abs.(kT * amhc.A .* pref .* erfc.(arg))
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

"""
    MarcusHushChidsey(A, λ, average_dos)
    MarcusHushChidsey(λ, average_dos)
    MarcusHushChidsey(λ)

Computes Marcus-Hush-Chidsey kinetics: 10.1126/science.251.4996.919

Note that for "typical" reorganization energy values (in the vicinity of 10*kT at typical temperatures, i.e. a few tenths of an eV), `AsymptoticMarcusHushChidsey` is comparably accurate to and much faster to evaluate than this model. 

If either the prefactor or the average dos are omitted, their values are assumed to be 1. Note that strictly speaking, `average_dos` and the prefactor `A` are redundant. They are both included primarily to facilitate comparisons with similarly parametrized Marcus-like models such as `MarcusHushChidseyDOS`.
"""
struct MarcusHushChidsey <: IntegralModel
    A::Real
    λ::Real
    average_dos::Real
end

# default prefactor is 1
MarcusHushChidsey(λ, avg_dos) = MarcusHushChidsey(1.0, λ, avg_dos)
# assume prefactor = 1 and avg_dos = 1
MarcusHushChidsey(λ) = MarcusHushChidsey(1.0, λ, 1.0)
# convert more detailed DOS information to just pull out average
MarcusHushChidsey(A, λ, dd::DOSData) = MarcusHushChidsey(A, λ, dd.average_value)
MarcusHushChidsey(A, λ, dos_file::Union{Matrix,String}; kwargs...) =
    MarcusHushChidsey(A, λ, DOSData(dos_file; kwargs...))

# TODO: Check that both this and +DOS versions still match original paper
function integrand(mhc::MarcusHushChidsey, V_dl, ox::Bool; kT::Real = 0.026)
    function marcus_term(E)
        local exp_arg
        if ox
            exp_arg = -((E .- mhc.λ .+ V_dl) .^ 2) ./ (4 * mhc.λ * kT)
        else
            exp_arg = -((E .- mhc.λ .- V_dl) .^2) ./ (4 * mhc.λ * kT)
        end
        exp.(exp_arg)
    end
    f(E::Vector) = hcat((mhc.A * mhc.average_dos) .* marcus_term.(E) .* fermi_dirac.(E; kT = kT)...)' # will return a matrix of both E and V_dl are vectors, first index will be energies and second voltages
    f(E::Real) = (mhc.A * mhc.average_dos) .* marcus_term.(E) .* fermi_dirac.(E; kT = kT)
end



"""
    MarcusHushChidseyDOS(A=1.0, λ, dos)
    MarcusHushChidseyDOS(A=1.0, λ, dos_file)

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

MarcusHushChidseyDOS(A, λ, dos_file::Union{Matrix,String}; kwargs...) =
    MarcusHushChidseyDOS(A, λ, DOSData(dos_file; kwargs...))

MarcusHushChidseyDOS(λ, dos_file::Union{Matrix,String}; kwargs...) = MarcusHushChidseyDOS(1.0, λ, DOSData(dos_file; kwargs...))

function integrand(
    mhcd::MarcusHushChidseyDOS,
    V_dl,
    ox::Bool;
    kT::Real = 0.026,
    V_q = 0.0,
)
    function marcus_term(E)
        local exp_arg
        if ox
            exp_arg = -(( (mhcd.λ .- V_dl) .+ E) .^ 2) ./ (4 * mhcd.λ * kT)
        else
            exp_arg = -(( (mhcd.λ .+ V_dl) .- E) .^ 2) ./ (4 * mhcd.λ * kT)
        end
        exp.(exp_arg)
    end
    fd(E) = ox ? 1 .- fermi_dirac(E; kT = kT) : fermi_dirac(E; kT = kT)
    # TODO: add two dispatches as with MHC above
    E -> mhcd.A .* mhcd.dos.interp_func.(E .+ V_q) .* fd(E) .* marcus_term(E)
end
