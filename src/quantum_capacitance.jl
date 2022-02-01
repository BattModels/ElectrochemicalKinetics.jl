using Interpolations
using QuadGK

function calculate_Vdl_interp(dos_f, Vq_min, Vq_max, C_dl)
    Vq_range = range(Vq_min, Vq_max, step=0.001)
    E_min = dos_f.ranges[1][1]
    E_max = dos_f.ranges[1][end]
    CQ = [compute_cq(E_min, E_max, V, dos_f)*1.602*10 for V in Vq_range]
    Vdl_data = []
    Vappl_data = []
    for i = 1:length(Vq_range)
        V_app_i =  Vq_range[i]*((CQ[i]/C_dl) + 1)
        V_dl_i = V_app_i - Vq_range[i]
        push!(Vdl_data, V_dl_i)
        push!(Vappl_data, V_app_i)
    end
    
    ## sort data and do V_dl interpolation
    inds = sortperm(Vappl_data)
    v_interp = LinearInterpolation(Vappl_data[inds], Vdl_data[inds])
    return v_interp
end


function QC_integrand(E, eVq, dos_func; kT=.026)
    dos_func(E) * sech((E-eVq)/(2*kT))^2/(4*kT)
end

function compute_cq(E_min, E_max, eVq, dos_func; kT=.026)
    fcn = E->QC_integrand(E, eVq, dos_func; kT=kT)
    quadgk(fcn, E_min, E_max)[1]
end
