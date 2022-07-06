"""
    ButlerVolmer(A, α)
    ButlerVolmer(A)

Computes Butler-Volmer kinetics. 

If initialized with one argument, assumes symmetric electron transfer (α=0.5) and sets this to be the prefactor A. Note that this prefactor implicitly contains information about equilibrium activation energies, as well as geometric information.
"""
struct ButlerVolmer{Ta, Talpha} <: NonIntegralModel
    A::Ta
    α::Talpha
end

function ButlerVolmer(A::T, α::S) where {T,S}
    @assert all(α .<= 1.0 .&& α .>=0.0) "Electron transfer coefficient must be in [0,1]"
    ButlerVolmer{T,S}(A, α)
end

# default to unit prefactor and symmetric response
ButlerVolmer() = ButlerVolmer(1.0, 0.5)
ButlerVolmer(α) = ButlerVolmer(1.0, α)

bv_f((V, ps), ::Val{true}; kT = 0.026) = ps[1] .* exp.((ps[2] .* V) / kT)
bv_f((V, ps), ::Val{false}; kT = 0.026) = ps[1] .* exp.(-((1 .- ps[2]) .* V) / kT)

function (bv::ButlerVolmer)(V, ::Val{true}; kT = convert(eltype(V), 0.026))
    bv.A .* exp.((bv.α .* V) ./ kT)
end

function (bv::ButlerVolmer)(V, ::Val{false}; kT = convert(eltype(V), 0.026))
    bv.A .* exp.(.-((one(eltype(V)) .- bv.α) .* V) ./ kT)
end

(bv::ButlerVolmer)(V; kT = convert(eltype(V), 0.026)) = abs.(bv(V, Val(true); kT) .- bv(V, Val(false); kT))

rate_f(::ButlerVolmer) = bv_f
