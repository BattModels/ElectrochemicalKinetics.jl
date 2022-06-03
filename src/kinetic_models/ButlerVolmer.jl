"""
    ButlerVolmer(A, α)
    ButlerVolmer(A)

Computes Butler-Volmer kinetics. 

If initialized with one argument, assumes symmetric electron transfer (α=0.5) and sets this to be the prefactor A. Note that this prefactor implicitly contains information about equilibrium activation energies, as well as geometric information.
"""
struct ButlerVolmer <: NonIntegralModel
    A
    α
end

# default to unit prefactor and symmetric response
ButlerVolmer() = ButlerVolmer(1.0, 0.5)
ButlerVolmer(A) = ButlerVolmer(A, 0.5)

bv_f((ps, V), ::Val{true}; kT = 0.026) = ps[1] .* exp.((ps[2] .* V) / kT)
bv_f((ps, V), ::Val{false}; kT = 0.026) = ps[1] .* exp.(-((1 .- ps[2]) .* V) / kT)

rate_f(::ButlerVolmer) = bv_f
