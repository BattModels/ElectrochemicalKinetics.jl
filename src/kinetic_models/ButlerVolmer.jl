"""
    ButlerVolmer(A, α)
    ButlerVolmer(A)

Computes Butler-Volmer kinetics. 

If initialized with one argument, assumes symmetric electron transfer (α=0.5) and sets this to be the prefactor A. Note that this prefactor implicitly contains information about equilibrium activation energies, as well as geometric information.
"""
struct ButlerVolmer{T} <: NonIntegralModel{T}
    A::T
    α::T
    function ButlerVolmer(A, α)
        @assert all(α .<= 1.0 .&& α .>=0.0) "Electron transfer coefficient must be in [0,1]"
        ps = consistent_params(Float32.(A), Float32.(α))
        new{typeof(ps[1])}(ps...)
    end
end

# default to unit prefactor and symmetric response
ButlerVolmer() = ButlerVolmer(1.0, 0.5)
ButlerVolmer(A) = ButlerVolmer(A, 0.5)

function rate_constant(V_app, bv::ButlerVolmer, ::Val{true}; kT::Real = 0.026)
    exp_arg = (bv.α .* V_app) / kT
    bv.A .* exp.(exp_arg)
end

function rate_constant(V_app, bv::ButlerVolmer, ::Val{false}; kT::Real = 0.026)
    exp_arg = -((1 .- bv.α) .* V_app) / kT
    bv.A .* exp.(exp_arg)
end
