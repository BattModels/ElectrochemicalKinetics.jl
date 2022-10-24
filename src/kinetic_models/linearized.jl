"""
    LinearizedKineticModel(A, T_f)

A simple linear kinetic model, used to make overpotential fitting more efficient near the origin. A linearized version of a "real" model can be constructed as `LinearizedKineticModel(m::KineticModel)`. The `T_f` field contains a function of a single variable (temperature) that encompasses the temperature dependence
"""

struct LinearizedKineticModel{T,KM} <: NonIntegralModel{T}
    A::T
    m_o::T
    m_r::T
    T_f::Function
    orig_model::KM
end

Base.show(io::IO, lkm::LinearizedKineticModel) = print(io, repr(typeof(lkm)))

LinearizedKineticModel(bv::ButlerVolmer) = LinearizedKineticModel(bv.A, bv.α, (bv.α-1), T->1, bv)

LinearizedKineticModel(m::Marcus) = LinearizedKineticModel(m.A, 0.5, -0.5, T->exp(-m.λ/(4*kB*T)), m)

rate_constant(V, lkm::LinearizedKineticModel, ::Val{true}; T=298) = lkm.A * lkm.T_f(T) .* (1 .+ lkm.m_o .* V ./ (kB*T))
rate_constant(V, lkm::LinearizedKineticModel, ::Val{false}; T=298) = lkm.A * lkm.T_f(T) .* (1 .+ lkm.m_r .* V ./ (kB*T))