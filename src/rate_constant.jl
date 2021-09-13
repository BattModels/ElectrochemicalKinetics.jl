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

Compute the rate constant k predicted by a given kinetic model at a applied voltage `V_app`. If a flag for reaction direction `ox` is supplied, `true` gives the oxidative and `false` the reductive direction, while omitting this flag will yield net reaction rate. If the model is an `IntegralModel`, integration bounds `E_min` and `E_max` must be supplied. Integration is done via GK quadrature.
"""
compute_k(V_app, model::KineticModel, ox::Bool; kwargs...) = model(V_app, ox; kwargs...)
compute_k(E_min, E_max, V_app, model::IntegralModel, ox::Bool; kwargs...) = quadgk(integrand(model, V_app, ox; kwargs...), E_min, E_max)[1] # return format is (value, error_bound)

compute_k(V_app, model::KineticModel; kwargs...) = model(V_app; kwargs...)
compute_k(E_min, E_max, V_app, model::IntegralModel; kwargs...) = quadgk(integrand(model, V_app; kwargs...), E_min, E_max)[1]

"""
    compute_k_cq(E_min, E_max, V_app, model::MarcusHushChidseyDOS, ox::Bool; kwargs...)
    compute_k_cq(E_min, E_max, V_app, model::MarcusHushChidseyDOS; kwargs...)

Compute the rate constant k predicted by a `MarcusHushChidseyDOS` model at a applied voltage `V_app`, including the effects of quantum capacitance. If a flag for reaction direction `ox` is supplied, `true` gives the oxidative and `false` the reductive direction, while omitting this flag will yield net reaction rate.
"""

function compute_k_cq(E_min, E_max, V_app, model::MarcusHushChidseyDOS, ox::Bool, C_dl=10.0, Vq_min=-0.5, Vq_max=0.5; kT = .026)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    V_dl = V_dl_interp(V_app)
    V_q = V_app - V_dl
    quadgk(integrand(model, V_dl, ox; kT=kT, V_q=V_q), E_min, E_max)[1]
end

function compute_k_cq(E_min, E_max, V_app, model::MarcusHushChidseyDOS, C_dl=10.0, Vq_min=-0.5, Vq_max=0.5; kT = .026)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    V_dl = V_dl_interp(V_app)
    V_q = V_app - V_dl
    quadgk(integrand(model, V_dl; kT=kT, V_q=V_q), E_min, E_max)[1]
end