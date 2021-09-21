using Statistics
using Interpolations
using QuadGK
using DelimitedFiles
using BlackBoxOptim
using Dierckx

# TODO: think about reaction direction and quantum cap dispatches; there's probably a cleaner way to handle them...

"""
    compute_k(V_app, model::KineticModel, ox::Bool; kwargs...)
    compute_k(V_app, model::KineticModel; kwargs...)
    compute_k(E_min, E_max, V_app, model::IntegralModel, ox::Bool; kwargs...)
    compute_k(E_min, E_max, V_app, model::IntegralModel; kwargs...)
    compute_k(E_min, E_max, V_app, model::MarcusHushChidseyDOS, calc_cq::Bool=false; C_dl = 10.0, Vq_min = -0.5, Vq_max = 0.5, kwargs...)

Compute the rate constant k predicted by a given kinetic model at a applied voltage `V_app`. If a flag for reaction direction `ox` is supplied, `true` gives the oxidative and `false` the reductive direction, while omitting this flag will yield net reaction rate. If the model is an `IntegralModel`, integration bounds `E_min` and `E_max` must be supplied. Integration is done via GK quadrature.

If calc_cq flag is passed, optionally compute voltage shifts due to quantum capacitance.
"""
compute_k(V_app, model::KineticModel, ox::Bool; kwargs...) = model(V_app, ox; kwargs...)
compute_k(E_min, E_max, V_app, model::IntegralModel, ox::Bool; kwargs...) =
    quadgk(integrand(model, V_app, ox; kwargs...), E_min, E_max)[1] # return format is (value, error_bound)

compute_k(V_app, model::KineticModel; kwargs...) = model(V_app; kwargs...)
compute_k(E_min, E_max, V_app, model::IntegralModel; kwargs...) =
    quadgk(integrand(model, V_app; kwargs...), E_min, E_max)[1]

function compute_k(E_min, E_max, V_app, model::MarcusHushChidseyDOS, ox::Bool; calc_cq::Bool=false, C_dl = 10.0, Vq_min = -0.5, Vq_max = 0.5, kwargs...)
    if calc_cq
        compute_k_cq(E_min, E_max, V_app, model, ox; C_dl=C_dl, Vq_min=Vq_min, Vq_max=Vq_max, kwargs...)
    else
        compute_k(E_min, E_max, V_app, model, ox; kwargs...)
    end
end

function compute_k(E_min, E_max, V_app, model::MarcusHushChidseyDOS; calc_cq::Bool=false, C_dl = 10.0, Vq_min = -0.5, Vq_max = 0.5, kwargs...)
    if calc_cq
        compute_k_cq(E_min, E_max, V_app, model, ox; C_dl=C_dl, Vq_min=VQ_min, Vq_max=Vq_max, kwargs...)
    else
        compute_k(E_min, E_max, V_app, model; kwargs...)
    end
end

"""
    compute_k_cq(E_min, E_max, V_app, model::MarcusHushChidseyDOS, ox::Bool; C_dl=10.0, Vq_min=-0.5, Vq_max=0.5, kwargs...)
    compute_k_cq(E_min, E_max, V_app, model::MarcusHushChidseyDOS; C_dl=10.0, Vq_min=-0.5, Vq_max=0.5, kwargs...)

Compute the rate constant k predicted by a `MarcusHushChidseyDOS` model at a applied voltage `V_app`, including the effects of quantum capacitance. If a flag for reaction direction `ox` is supplied, `true` gives the oxidative and `false` the reductive direction, while omitting this flag will yield net reaction rate.
"""
function compute_k_cq(
    E_min,
    E_max,
    V_app,
    model::MarcusHushChidseyDOS,
    ox::Bool;
    C_dl = 10.0,
    Vq_min = -0.5,
    Vq_max = 0.5,
    kT = 0.026,
)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    V_dl = V_dl_interp(V_app)
    V_q = V_app - V_dl
    quadgk(integrand(model, V_dl, ox; kT = kT, V_q = V_q), E_min, E_max)[1]
end

function compute_k_cq(
    E_min,
    E_max,
    V_app,
    model::MarcusHushChidseyDOS;
    C_dl = 10.0,
    Vq_min = -0.5,
    Vq_max = 0.5,
    kT = 0.026,
)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    V_dl = V_dl_interp(V_app)
    V_q = V_app - V_dl
    quadgk(integrand(model, V_dl; kT = kT, V_q = V_q), E_min, E_max)[1]
end
