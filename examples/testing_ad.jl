using BenchmarkTools
using Zygote
using ElectrochemicalKinetics

# THE CURRENT SETUP

# benchmark just the forward evaluation (~400 ns on my laptop):
@benchmark bv($0.1, $Val(true))

# and just the analytical gradient (also about 400 ns)
dbv(V, bv::ButlerVolmer, ::Val{true}, kT=.026) = bv.α/kT * bv(V)
@benchmark dbv($0.1, $bv, $Val(true))

# now for gradients
# you can try this next line with and without the adjoint in lines 63-67 of kinetic_models.jl commented out, it goes from ~5 µs to ~4 µs and from ~3.4 KiB allocations to ~1.7
@benchmark gradient(V->bv(V, $Val(true)), $0.1)

# ALTERNATIVE SETUP
# and here's a test of a version where everything is functions as opposed to this callable struct setup, which seems to be a lot faster...
bvf(V, ox=true; α=0.5, A=20, kT=0.026) = A * exp(α * V / kT)
@benchmark bvf($0.1) # about 10 ns, no allocations
dbv(V, ox=true; α=0.5, A=20, kT=0.026) = A * α / kT * exp(α * V / kT) 
@benchmark dbv($0.1) # similar
@adjoint bvf(V, ox; kwargs...) = bvf(V, ox; kwargs...), x -> (x * dbv(V, ox; kwargs...), nothing)

# ~16 ns, no allocations (so actually less than the sum of the two above, rather than 10x like before)
@benchmark gradient(V->bvf(V), $0.1)

# REALISTIC VERSION with compute_k
# compute_k(V_app, model::KineticModel, args...; kT = 0.026) = model(V_app, args...; kT = kT)
# for now only for ox=true case...
compute_k(V_app, bv::ButlerVolmer; kT = 0.026) = bv.A * exp(bv.α * V_app / kT)
@benchmark compute_k($0.01, $bv) # about 10 ns, no allocations

@adjoint function compute_k(V_app, bv::ButlerVolmer; kT = 0.026)
    k = compute_k(V_app, bv; kT=kT)
    k, x -> (x * k * bv.α / kT, nothing)
end

@benchmark gradient(V->compute_k(V, $bv), $0.1) # about 9 ns