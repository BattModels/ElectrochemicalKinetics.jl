using Statistics
using Interpolations
using QuadGK
using DelimitedFiles
using BlackBoxOptim
using Dierckx



"""
    compute_k(E_min, E_max, λ, V_dl, model)
"""
function compute_k(E_min, E_max, model::KineticModel, direction; kT=.026, C_q_calc=false, C_dl=10.0, Vq_min=-0.5, Vq_max=0.5)
    @assert direction in ["red", "ox", "net"] "The supported reaction directions are red, ox, and net"

    
end

"""
    compute_k_MHC(E_min, E_max, λ, V_dl, average_dos, ox=true; kT=.026)

Compute the Marcus-Hush-Chidsey rate constant via GK quadrature by integrating the result of [`MHC_integrand`](@ref)(E, λ, V_dl, average_dos, ox; kT) from `E`=`E_min` to `E_max`.
"""
function compute_k_MHC(E_min, E_max, λ, V_dl, average_dos, ox; kT=.026)
    fcn = E -> MHC_integrand(E, λ, V_dl, average_dos, ox; kT=.026)
    quadgk(fcn, E_min, E_max)[1] # return format is (value, error_bound)
end

"""
    compute_k_MHC_net(E_min, E_max, λ, V_dl, average_dos; kT=.026)

Compute the net Marcus-Hush-Chidsey rate constant via GK quadrature by integrating the result of [`MHC_integrand_net`](@ref)(E, λ, V_dl, average_dos; kT) from `E`=`E_min` to `E_max`.
"""
function compute_k_MHC_net(E_min, E_max, λ, V_dl, average_dos; kT=.026)
    fcn = E -> MHC_integrand_net(E, λ, V_dl, average_dos; kT=.026)
    quadgk(fcn, E_min, E_max)[1] # return format is (value, error_bound)
end


function compute_k_MHC_DOS(E_min, E_max, λ, V_app, dos_func; kT=.026, C_q_calc=false, C_dl=10.0, Vq_min=-0.5, Vq_max=0.5)
    if C_q_calc
        ## calculate V_dl from V_app by calculating CQ
        V_dl_interp = calculate_Vdl_interp(dos_func, Vq_min, Vq_max, C_dl)
        V_dl = V_dl_interp(V_app)
        V_q = V_app - V_dl
        fcn = E->MHC_DOS_integrand_net(E, λ, V_dl, dos_func; kT=kT, V_q=V_q)
    else
        fcn = E->MHC_DOS_integrand_net(E, λ, V_app, dos_func; kT=kT, V_q=0)
    end
    quadgk(fcn, E_min, E_max)[1]
end


function compute_k_MHC_DOS_redox(E_min, E_max, λ, V_app, dos_func; ox=true, kT=.026, C_q_calc=false, C_dl=10.0, Vq_min=-0.5, Vq_max=0.5)
    if C_q_calc
        V_dl_interp = calculate_Vdl_interp(dos_func, Vq_min, Vq_max, C_dl)
        V_dl = V_dl_interp(V_app)
        V_q = V_app - V_dl
        fcn = E->MHC_DOS_integrand(E, λ, V_dl, dos_func; ox=ox, kT=kT, V_q=V_q)
    else
        fcn = E->MHC_DOS_integrand(E, λ, V_dl, dos_func; ox=ox, kT=kT, V_q=0)
    end
    quadgk(fcn, E_min, E_max)[1]
end


