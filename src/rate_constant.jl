using Statistics
using Interpolations
# using QuadGK
using DelimitedFiles
using Dierckx

# TODO: think about reaction direction and quantum cap dispatches; there's probably a cleaner way to handle them...

"""
    compute_k(V_app, model::KineticModel, ox::Bool; kwargs...)
    compute_k(V_app, model::KineticModel; kwargs...)
    compute_k(E_min, E_max, V_app, model::IntegralModel, ox::Bool; kwargs...)
    compute_k(E_min, E_max, V_app, model::IntegralModel; kwargs...)
    compute_k(E_min, E_max, V_app, model::MarcusHushChidseyDOS, calc_cq::Bool=false; C_dl = 10.0, Vq_min = -0.5, Vq_max = 0.5, kwargs...)

Compute the rate constant k predicted by a given kinetic model at a applied voltage `V_app`. If a flag for reaction direction `ox` is supplied, `true` gives the oxidative and `false` the reductive direction, while omitting this flag will yield net reaction rate (absolute value thereof).

If the model is an `IntegralModel`, integration bounds `E_min` and `E_max` must be supplied. Integration is done via GK quadrature.

If calc_cq flag is passed, optionally compute voltage shifts due to quantum capacitance.
"""
compute_k(V_app, model::KineticModel, ox::Bool; kT = 0.026) = model(V_app, ox; kT = kT)

compute_k(
    V_app,
    model::IntegralModel,
    ox::Bool;
    kT = 0.026,
    E_min = -100 * kT,
    E_max = 100 * kT,
) = begin
  n, w = scale(E_min, E_max)
  f = integrand(model, V_app, ox; kT = kT)
  sum(w .* f(n))
  # quadgk(integrand(model, V_app, ox; kT = kT), E_min, E_max)[1] # return format is (value, error_bound)
end

compute_k(V_app, model::KineticModel; kT = 0.026) = model(V_app; kT = kT)
compute_k(
    V_app,
    model::IntegralModel;
    kT = 0.026,
    E_min = -100 * kT,
    E_max = 100 * kT,
    kwargs...,
) = begin
  n, w = scale(E_min, E_max)
  f = integrand(model, V_app; kT = kT)
  sum(w .* f.(n))
  # quadgk(integrand(model, V_app; kT = kT), E_min, E_max)[1]
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
        # quadgk(integrand(model, V_app, ox; kT = kT), E_min, E_max)[1]
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
        # quadgk(integrand(model, V_app; kT = kT), E_min, E_max)[1]
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
    # quadgk(integrand(model, V_dl, ox; kT = kT, V_q = V_q), E_min, E_max)[1]
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
    # quadgk(integrand(model, V_dl; kT = kT, V_q = V_q), E_min, E_max)[1]
end
