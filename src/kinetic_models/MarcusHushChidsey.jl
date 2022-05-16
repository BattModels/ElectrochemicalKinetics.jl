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