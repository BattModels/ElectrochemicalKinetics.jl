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