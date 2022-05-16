include("../utils/quantum_capacitance.jl")

"""
    MarcusHushChidseyDOS(A=1.0, λ, dos)
    MarcusHushChidseyDOS(A=1.0, λ, dos_file)

Computes Marcus-Hush-Chidsey + DOS kinetics as described in Kurchin and Viswanathan: 10.1063/5.0023611 
"""
struct MarcusHushChidseyDOS <: IntegralModel
    A::Float64
    λ::Float64
    dos::DOSData
end

function Base.show(io::IO, mhcd::MarcusHushChidseyDOS)
    s = repr(typeof(mhcd)) * "("
    s *= "A=$(round(mhcd.A, sigdigits=3)), λ=$(round(mhcd.λ, sigdigits=3)))"
    print(io, s)
end

# default prefactor is 1
MarcusHushChidseyDOS(λ, dd::DOSData) = MarcusHushChidseyDOS(1.0, λ, dd)

MarcusHushChidseyDOS(A, λ, dos_file::Union{Matrix,String}; kwargs...) =
    MarcusHushChidseyDOS(A, λ, DOSData(dos_file; kwargs...))

MarcusHushChidseyDOS(λ, dos_file::Union{Matrix,String}; kwargs...) = MarcusHushChidseyDOS(1.0, λ, DOSData(dos_file; kwargs...))

function integrand(
    mhcd::MarcusHushChidseyDOS,
    V_dl,
    ox::Bool;
    kT::Real = 0.026,
    V_q = 0.0,
)
    function marcus_term(E)
        local exp_arg
        if ox
            exp_arg = -(( (mhcd.λ .- V_dl) .+ E) .^ 2) ./ (4 * mhcd.λ * kT)
        else
            exp_arg = -(( (mhcd.λ .+ V_dl) .- E) .^ 2) ./ (4 * mhcd.λ * kT)
        end
        exp.(exp_arg)
    end
    fd(E) = ox ? 1 .- fermi_dirac(E; kT = kT) : fermi_dirac(E; kT = kT)
    # TODO: add two dispatches as with MHC above
    E -> mhcd.A .* mhcd.dos.interp_func.(E .+ V_q) .* fd(E) .* marcus_term(E)
end

function compute_k(
    V_app,
    model::MarcusHushChidseyDOS,
    ox::Bool;
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
        compute_k_cq(
            V_app,
            model,
            ox;
            C_dl = C_dl,
            Vq_min = Vq_min,
            Vq_max = Vq_max,
            kT = kT,
            E_min = min(E_min, E_min-Vq_min),
            E_max = max(E_max, E_max-Vq_max),
        )
    else
        n, w = scale(E_min, E_max)
        f = integrand(model, V_app, ox; kT = kT)
        sum(w .* f.(n))
    end
end

function compute_k(
    V_app,
    model::MarcusHushChidseyDOS;
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
        compute_k_cq(
            V_app,
            model;
            kT = kT,
            E_min = max(E_min, E_min-Vq_min),
            E_max = min(E_max, E_max-Vq_max),
            C_dl = C_dl,
            Vq_min = Vq_min,
            Vq_max = Vq_max,
        )
    else
        # invoke(
        #     compute_k,
        #     Tuple{typeof(V_app),IntegralModel},
        #     V_app,
        #     model;
        #     kT = kT,
        #     E_min = E_min,
        #     E_max = E_max,
        # )
        n, w = scale(E_min, E_max)
        f = integrand(model, V_app; kT = kT)
        sum(w .* f.(n))
    end
end

"""
    compute_k_cq(E_min, E_max, V_app, model::MarcusHushChidseyDOS, ox::Bool; C_dl=10.0, Vq_min=-0.5, Vq_max=0.5, kwargs...)
    compute_k_cq(E_min, E_max, V_app, model::MarcusHushChidseyDOS; C_dl=10.0, Vq_min=-0.5, Vq_max=0.5, kwargs...)

Compute the rate constant k predicted by a `MarcusHushChidseyDOS` model at a applied voltage `V_app`, including the effects of quantum capacitance. If a flag for reaction direction `ox` is supplied, `true` gives the oxidative and `false` the reductive direction, while omitting this flag will yield net reaction rate.
"""
function compute_k_cq(
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

function compute_k_cq(
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
