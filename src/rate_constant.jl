using Statistics
using Interpolations
using QuadGK
using DelimitedFiles
using BlackBoxOptim
using Dierckx

# TODO: think about reaction direction and quantum cap dispatches; there's probably a cleaner way to handle them...

# for one direction
compute_k(E_min, E_max, V_app, model::KineticModel, ox::Bool; args...) = quadgk(integrand(model, V_app, ox; args...), E_min, E_max)[1] # return format is (value, error_bound)

# for net
compute_k(E_min, E_max, V_app, model::KineticModel; args...) = quadgk(integrand(model, V_app; args...), E_min, E_max)[1]

function compute_k_cq(E_min, E_max, V_app, model::MarcusHushChidseyDOS, ox::Bool, C_dl=10.0, Vq_min=-0.5, Vq_max=0.5; kT = .026)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    V_dl = V_dl_interp(V_app)
    V_q = V_app - V_dl
    quadgk(integrand(model, V_dl, ox; kT=kT, V_q=V_q), E_min, E_max)[1]
end

# TODO: check that this gives the same as computing them separately and subtracting in all cases (including with V_q)
function compute_k_cq(E_min, E_max, V_app, model::MarcusHushChidseyDOS, C_dl=10.0, Vq_min=-0.5, Vq_max=0.5; kT = .026)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    V_dl = V_dl_interp(V_app)
    V_q = V_app - V_dl
    quadgk(integrand(model, V_dl; kT=kT, V_q=V_q), E_min, E_max)[1]
end