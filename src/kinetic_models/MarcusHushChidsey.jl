"""
    MarcusHushChidsey(A, λ, average_dos)
    MarcusHushChidsey(λ, average_dos)
    MarcusHushChidsey(λ)

Computes Marcus-Hush-Chidsey kinetics: 10.1126/science.251.4996.919

Note that for "typical" reorganization energy values (in the vicinity of 10*kT at typical temperatures, i.e. a few tenths of an eV), `AsymptoticMarcusHushChidsey` is comparably accurate to and much faster to evaluate than this model. 

If either the prefactor or both the prefactor and the average dos are omitted, their values are assumed to be 1. Note that strictly speaking, `average_dos` and the prefactor `A` are redundant. They are both included primarily to facilitate comparisons with similarly parametrized Marcus-like models such as `MarcusHushChidseyDOS`.
"""
struct MarcusHushChidsey <: IntegralModel
    A
    λ
    average_dos
end

# default prefactor is 1
MarcusHushChidsey(λ, avg_dos) = MarcusHushChidsey(1.0, λ, avg_dos)
# assume prefactor = 1 and avg_dos = 1
MarcusHushChidsey(λ) = MarcusHushChidsey(1.0, λ, 1.0)
# convert more detailed DOS information to just pull out average
MarcusHushChidsey(A, λ, dd::DOSData) = MarcusHushChidsey(A, λ, Ref(dd.average_value))
MarcusHushChidsey(A, λ, dos_file::Union{Matrix,String}; kwargs...) =
    MarcusHushChidsey(A, λ, DOSData(dos_file; kwargs...))

# TODO: Check that both this and +DOS versions still match original paper
function integrand(mhc::MarcusHushChidsey, V_dl, ox::Val; T::Real=298)
    kT = kB * T
    mhc_f((E, V, ps), ::Val{true}) = ps[1] * exp.(-((E .- ps[2] .+ V) .^ 2) ./ (4 .* ps[2] .* kT)) * fermi_dirac(E; T=T)
    mhc_f((E, V, ps), ::Val{false}) = ps[1] * exp.(-((E .- ps[2] .- V) .^2) ./ (4 .* ps[2] .* kT)) * fermi_dirac(E; T=T)
    pref = mhc.A .* mhc.average_dos
    # TODO: decide if this ordering makes sense (possibly flip model params and V?)
    f(E) = mhc_f.(Iterators.product(E, V_dl, iterate_props(mhc)), Ref(ox))
end
