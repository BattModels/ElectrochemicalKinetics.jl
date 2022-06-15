"""
    Marcus(A, λ)
    Marcus(λ)

Computes Marcus kinetics.

If initialized with one argument, assumes this to be the reorganization energy λ and sets the prefactor A to 1.0.
"""
struct Marcus <: NonIntegralModel
    A
    λ
end

# default prefactor is 1
Marcus(λ) = Marcus(1.0, λ)

m_f((V, ps), ::Val{true}; kT = 0.026) = ps[1] .* exp.(-(ps[2].+ V).^2 ./ (4 .* ps[2] .* kT))
m_f((V, ps), ::Val{false}; kT = 0.026) = ps[1] .* exp.(-(ps[2].- V).^2 ./ (4 .* ps[2] .* kT))

rate_f(::Marcus) = m_f
