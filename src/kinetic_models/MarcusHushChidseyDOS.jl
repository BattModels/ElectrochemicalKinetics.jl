include("../utils/quantum_capacitance.jl")

"""
    MarcusHushChidseyDOS(A=1.0, λ, dos)
    MarcusHushChidseyDOS(A=1.0, λ, dos_file)

Computes Marcus-Hush-Chidsey + DOS kinetics as described in Kurchin and Viswanathan: 10.1063/5.0023611

NB: At the moment, this will allow for vector `A` and `λ` parameters, but will presume that all models correspond to the same DOS.
"""
struct MarcusHushChidseyDOS <: IntegralModel
    A
    λ
    dos::DOSData
end

function Base.show(io::IO, mhcd::MarcusHushChidseyDOS)
    s = repr(typeof(mhcd)) * "("
    s *= "A=$(round.(mhcd.A, sigdigits=3)), λ=$(round(mhcd.λ, sigdigits=3)))"
    print(io, s)
end

# default prefactor is 1
MarcusHushChidseyDOS(λ, dd::DOSData) = MarcusHushChidseyDOS(1.0, λ, dd)

MarcusHushChidseyDOS(A, λ, dos_file::Union{Matrix,String}; kwargs...) =
    MarcusHushChidseyDOS(A, λ, DOSData(dos_file; kwargs...))

MarcusHushChidseyDOS(λ, dos_file::Union{Matrix,String}; kwargs...) = MarcusHushChidseyDOS(1.0, λ, DOSData(dos_file; kwargs...))

# TODO: make one version of this for the whole package, currently there are some sign convention differences between this and the MHC version 
marcus_term_ox(V, E, λ, kT) = exp.(-(( (λ .- V) .+ E) .^ 2) ./ (4 * λ * kT))
marcus_term_red(V, E, λ, kT) = exp.(-(( (λ .+ V) .- E) .^ 2) ./ (4 * λ * kT))

# Fermi-Dirac function or the complement thereof
fd_ox(E, kT) = 1 .- fermi_dirac(E; kT = kT)
fd_red(E, kT) = fermi_dirac(E; kT = kT)

function integrand(
    mhcd::MarcusHushChidseyDOS,
    V_dl,
    ox::Val;
    kT = 0.026,
    V_q = 0.0,
)
    mhcd_f((E, V, ps), ::Val{true}) = ps[1] * mhcd.dos.interp_func(E + V_q) * fd_ox(E, kT) * marcus_term_ox(V, E, ps[2], kT)
    mhcd_f((E, V, ps), ::Val{false}) = ps[1] * mhcd.dos.interp_func(E + V_q) * fd_red(E, kT) * marcus_term_red(V, E, ps[2], kT)

    f(E) = mhcd_f.(Iterators.product(E, V_dl, iterate_props(mhcd, exclude_props = [:dos])), Ref(ox))
end

function rate_constant(
    V_app,
    model::MarcusHushChidseyDOS,
    args...;
    kT = 0.026,
    E_min = model.dos.E_min,
    E_max = model.dos.E_max,
    calc_cq::Bool = false,
    C_dl = 10.0,
    Vq_min = -0.5,
    Vq_max = 0.5,
    kwargs...,
)
    if calc_cq
        res = rate_constant_cq(
            V_app,
            model,
            args...;
            C_dl = C_dl,
            Vq_min = Vq_min,
            Vq_max = Vq_max,
            kT = kT,
            E_min = min(E_min, E_min-Vq_min),
            E_max = max(E_max, E_max-Vq_max),
        )
    else
        n, w = scale(E_min, E_max)
        f = integrand(model, V_app, args...; kT = kT)
        res = sum(w .* f.(n))
    end

    if size(res) == (1,)
        return res[1]
    else
        return res
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
    kT = 0.026,
    E_min = model.dos.E_min,
    E_max = model.dos.E_max,
)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    V_dl = V_dl_interp(V_app)
    V_q = V_app - V_dl
    n, w = scale(E_min, E_max)
    f = integrand(model, V_dl, ox; kT = kT, V_q = V_q)
    sum(w .* f.(n))
end

# TODO: clean these up...but add tests first
function rate_constant_cq(
    V_app,
    model::MarcusHushChidseyDOS;
    C_dl = 10.0,
    Vq_min = -0.5,
    Vq_max = 0.5,
    kT = 0.026,
    E_min = model.dos.E_min,
    E_max = model.dos.E_max,
)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    V_dl = V_dl_interp(V_app)
    V_q = V_app - V_dl
    n, w = scale(E_min, E_max)
    f = integrand(model, V_dl; kT = kT, V_q = V_q)
    sum(w .* f.(n))
end
