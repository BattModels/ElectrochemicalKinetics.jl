"""
    Marcus(A, λ)
    Marcus(λ)

Computes Marcus kinetics.

If initialized with one argument, assumes this to be the reorganization energy λ and sets the prefactor A to 1.0.
"""
struct Marcus{T} <: NonIntegralModel{T}
    A::T
    λ::T
    function Marcus(A, λ)
        ps = consistent_params(Float64.(A), Float64.(λ))
        new{typeof(ps[1])}(ps...)
    end
end

# default prefactor is 1
Marcus(λ) = Marcus(1.0, λ)

function rate_constant(V_app, m::Marcus, ::Val{true}; kT = 0.026)
    exp_arg = -(m.λ .+ V_app).^2 ./ (4 .* m.λ .* kT)
    m.A .* exp.(exp_arg)
end

function rate_constant(V_app, m::Marcus, ::Val{false}; kT = 0.026)
    exp_arg = -(m.λ .- V_app).^2 ./ (4 .* m.λ .* kT)
    m.A .* exp.(exp_arg)
end