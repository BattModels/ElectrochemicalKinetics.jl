"""
    ButlerVolmer(A, α)
    ButlerVolmer(A)

Computes Butler-Volmer kinetics. 

If initialized with one argument, assumes symmetric electron transfer (α=0.5) and sets this to be the prefactor A. Note that this prefactor implicitly contains information about equilibrium activation energies, as well as geometric information.
"""
struct ButlerVolmer <: NonIntegralModel
    A
    α
    function ButlerVolmer(A, α)
        @assert all(α .<= 1.0 .&& α .>=0.0) "Electron transfer coefficient must be in [0,1]"
        new(A, α)
    end
end

# default to unit prefactor and symmetric response
ButlerVolmer() = ButlerVolmer(1.0, 0.5)
ButlerVolmer(α) = ButlerVolmer(1.0, α)

bv_f((V, ps), ::Val{true}; T=298) = ps[1] .* exp.((ps[2] .* V) / (kB*T))
bv_f((V, ps), ::Val{false}; T = 298) = ps[1] .* exp.(-((1 .- ps[2]) .* V) / (kB*T))

rate_f(::ButlerVolmer) = bv_f
