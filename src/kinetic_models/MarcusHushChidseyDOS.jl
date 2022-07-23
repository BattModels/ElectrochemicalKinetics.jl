include("../utils/quantum_capacitance.jl")

"""
    MarcusHushChidseyDOS(A=1.0, λ, dos)
    MarcusHushChidseyDOS(A=1.0, λ, dos_file)

Computes Marcus-Hush-Chidsey + DOS kinetics as described in Kurchin and Viswanathan: 10.1063/5.0023611

NB: At the moment, this will allow for vector `A` and `λ` parameters, but will presume that all models correspond to the same DOS.
"""
struct MarcusHushChidseyDOS{T} <: IntegralModel{T}
    A::T
    λ::T
    dos::DOSData
    function MarcusHushChidseyDOS(A, λ, dos)
        ps = consistent_params(Float32.(A), Float32.(λ))
        new{typeof(ps[1])}(ps..., dos)
    end
end

function Base.show(io::IO, mhcd::MarcusHushChidseyDOS)
    s = repr(typeof(mhcd)) * "("
    s *= "A=$(round.(mhcd.A, sigdigits=3)), λ=$(round.(mhcd.λ, sigdigits=3)))"
    print(io, s)
end

# default prefactor is 1
MarcusHushChidseyDOS(λ, dd::DOSData) = MarcusHushChidseyDOS(1.0, λ, dd)

MarcusHushChidseyDOS(A, λ, dos_file::Union{Matrix,String}; kwargs...) =
    MarcusHushChidseyDOS(A, λ, DOSData(dos_file; kwargs...))

MarcusHushChidseyDOS(λ, dos_file::Union{Matrix,String}; kwargs...) = MarcusHushChidseyDOS(1.0, λ, DOSData(dos_file; kwargs...))

# TODO: make one version of this for the whole package, currently there are some sign convention differences between this and the MHC version 
marcus_term(V, E, λ, ::Val{true}, T) = exp.(-(( (λ .- V) .+ E) .^ 2) ./ (4 * λ * kB * T))
marcus_term(V, E, λ, ::Val{false}, T) = exp.(-(( (λ .+ V) .- E) .^ 2) ./ (4 * λ * kB * T))

# Fermi-Dirac function or the complement thereof
fd(E, ::Val{true}, T) = 1 .- fermi_dirac(E; T = T)
fd(E, ::Val{false}, T) = fermi_dirac(E; T = T)

function integrand(
    mhcd::MarcusHushChidseyDOS,
    V_dl,
    ox::Val;
    T = 298,
    V_q = 0.0,
)
    E -> mhcd.A .* mhcd.dos.interp_func.(E .+ V_q) .* fd(E, ox, T) .* marcus_term(V_dl, E, mhcd.λ, ox, T)
end

function rate_constant(
    V_app,
    model::MarcusHushChidseyDOS,
    args...;
    T = 298,
    E_min = model.dos.E_min,
    E_max = model.dos.E_max,
    calc_cq::Bool = false,
    C_dl = 10.0,
    Vq_min = -0.5,
    Vq_max = 0.5,
    kwargs...,
)
    if calc_cq
        rate_constant_cq(
            V_app,
            model,
            args...;
            C_dl = C_dl,
            Vq_min = Vq_min,
            Vq_max = Vq_max,
            T = T,
            E_min = min(E_min, E_min-Vq_min),
            E_max = max(E_max, E_max-Vq_max),
        )
    else
        n, w = scale(E_min, E_max)
        f = integrand(model, V_app, args...; T = T)
        sum(w .* f.(n))
    end
end

"""
    rate_constant_cq(E_min, E_max, V_app, model::MarcusHushChidseyDOS, ox::Bool; C_dl=10.0, Vq_min=-0.5, Vq_max=0.5, kwargs...)
    rate_constant_cq(E_min, E_max, V_app, model::MarcusHushChidseyDOS; C_dl=10.0, Vq_min=-0.5, Vq_max=0.5, kwargs...)

Compute the rate constant k predicted by a `MarcusHushChidseyDOS` model at a applied voltage `V_app`, including the effects of quantum capacitance. If a flag for reaction direction `ox` is supplied, `true` gives the oxidative and `false` the reductive direction, while omitting this flag will yield net reaction rate.
"""
function rate_constant_cq(
    V_app,
    model::MarcusHushChidseyDOS,
    ox::Bool;
    C_dl = 10.0,
    Vq_min = -0.5,
    Vq_max = 0.5,
    T = 298,
    E_min = model.dos.E_min,
    E_max = model.dos.E_max,
)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    V_dl = V_dl_interp(V_app)
    V_q = V_app - V_dl
    n, w = scale(E_min, E_max)
    f = integrand(model, V_dl, ox; T = T, V_q = V_q)
    sum(w .* f.(n))
end

# TODO: clean these up...but add tests first
function rate_constant_cq(
    V_app,
    model::MarcusHushChidseyDOS;
    C_dl = 10.0,
    Vq_min = -0.5,
    Vq_max = 0.5,
    T = 298,
    E_min = model.dos.E_min,
    E_max = model.dos.E_max,
)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    V_dl = V_dl_interp(V_app)
    V_q = V_app - V_dl
    n, w = scale(E_min, E_max)
    f = integrand(model, V_dl; T = T, V_q = V_q)
    sum(w .* f.(n))
end
