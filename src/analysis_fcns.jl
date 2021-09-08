using Statistics
using Interpolations
using QuadGK
#using LaTeXStrings
using DelimitedFiles
using BlackBoxOptim
using Dierckx

"""
    fermi_dirac(E; kT=0.026)

Compute the value of the Fermi-Dirac distribution at energy `E` (relative to the Fermi energy) and thermal energy `kT`.
"""
fermi_dirac(E; kT=.026) = 1 / (1 + exp(E/kT))

"""
    get_dos(dos_file; Ef=0, cut_energy=false)

Retrieve density of states data from a file at path `dos_file`, assuming that the Fermi level is at energy `Ef`.

Returns a callable interpolated DOS, the average value of DOS, and the lower and upper bound energies.

# Notes
* `dos_file` is assumed to contain two colums, energy and DOS value
* If `cut_energy==true`, return data in a symmetric range of energies around `Ef` (this is sometimes useful for nice-looking plots). Note that the code presumes the cutoff for this comes from the upper bound. This would not be hard to fix, I just have not done so as yet.
"""
function get_dos(dos_file; Ef=0, cut_energy=false)
    dos_data = readdlm(dos_file, Float64, skipstart=1)
    # recenter so Ef=0
    dos_data[:,1] = dos_data[:,1] .- Ef
    E_max = dos_data[end,1]
    if cut_energy #TODO: fix this so that it works if magnitude of upper bound is smaller
        # for now, we'll pick a symmetric energy range about 0 (i.e. the new Ef)
        len_keep = sum(dos_data[:,1] .> -E_max)
        # cut off
        dos_data = dos_data[end-len_keep:end,:]
    end
    E_min = dos_data[1,1]
    average_dos = mean(dos_data[:,2]) # for whole structure
    # interpolate
    E_step = mean(dos_data[2:end,1].-dos_data[1:end-1,1])
    dos_interp = scale(interpolate(dos_data[:,2], BSpline(Linear())), range(E_min, E_max+0.0001, step=E_step))
    return dos_interp, average_dos, E_min, E_max
end

"""
    marcus_integrand(E, λ, V_dl, ox=true; kT=.026)

Compute the value of the Marcus integrand, given by:
``-(λ - V_{dl} + E)^2 / (4 λ kT)``
in the oxidative direction (when `ox==true`), and
``-(λ + V_{dl} - E)^2 / (4 λ kT)``
in the reductive direction (when `ox==false`)

# Arguments (all energies are in eV)
* `E`: energy coordinate at which to evaluate integrand
* `λ`: reorganization energy parameter
* `V_dl`: voltage at the interface (equal to applied voltage in absence of quantum capacitance, see [`calculate_Vdl_interp`](@ref))
* `ox`: boolean flag for oxidative vs. reductive direction of reaction
* `kT`: thermal energy, defaults to room temperature value
"""
function marcus_integrand(E, λ, V_dl, ox=true; kT=.026)
    if ox # oxidative direction
        arg = -(λ-V_dl+E)^2 / (4*λ*kT)
    else # reductive direction
        arg = -(λ+V_dl-E)^2 / (4*λ*kT)
    end
    exp(arg)
end

"""
    MHC_integrand(E, λ, V_dl, average_dos, ox=true; kT=.026)

Compute the value of the Marcus-Hush-Chidsey integrand, given by the value of the Marcus integrand (see [`marcus_integrand`](@ref)) multiplied by the appropriate Fermi-Dirac distribution (for electrons or holes).

# Arguments (all energies are in eV)
* `E`: energy coordinate at which to evaluate integrand
* `λ`: reorganization energy parameter
* `V_dl`: voltage at the interface (equal to applied voltage in absence of quantum capacitance, see [`calculate_Vdl_interp`](@ref))
* `average_dos`: average value of the density of states (this is merely a scaling factor, useful for direct comparisons to the results of the full MHC+DOS method)
* `ox`: boolean flag for oxidative vs. reductive direction of reaction
* `kT`: thermal energy, defaults to room temperature value
"""
function MHC_integrand(E, λ, V_dl, average_dos; ox=true, kT=.026)
    marcus = marcus_integrand(E, λ, V_dl, ox; kT=kT)
    fd = ox ? 1-fermi_dirac(E; kT=kT) : fermi_dirac(E; kT=kT)
    average_dos * marcus * fd
end

# TODO: rewrite this in terms of ox and red fcns
function MHC_integrand_net(E, λ, V_dl, average_dos; kT=.026)
    marcus_ox = marcus_integrand(E, λ, V_dl, true; kT=kT)
    marcus_red = marcus_integrand(E, λ, V_dl, false; kT=kT)
    fd_ox = 1-fermi_dirac(E; kT=kT) # holes
    fd_red = fermi_dirac(E; kT=kT) # electrons
    sign(V_dl) * average_dos * (marcus_ox*fd_ox - marcus_red*fd_red)
end


function compute_k_MHC(E_min, E_max, λ, V_dl, average_dos; kT=.026)
    fcn = E -> MHC_integrand_net(E, λ, V_dl, average_dos; kT=.026)
    quadgk(fcn, E_min, E_max)[1] # return format is (value, error_bound)
end

# compute with DOS by evaluating the function we already have for a unit DOS and multiplying
# by the interpolated value
function MHC_DOS_integrand(E, λ, V_dl, dos_func; kT=.026, V_q=0)
    dos_func(E+V_q) * MHC_integrand_net(E, λ, V_dl, 1; kT=kT)
end

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
    
    ## do V_dl interpolation
    v_interp = LinearInterpolation(Vappl_data, Vdl_data)
    return v_interp
end

# TODO: rewrite this in terms of redox fcns
function compute_k_MHC_DOS(E_min, E_max, λ, V_app, dos_func; kT=.026, C_q_calc=false, C_dl=10.0, Vq_min=-0.5, Vq_max=0.5)
    if C_q_calc
        ## calculate V_dl from V_app by calculating CQ
        V_dl_interp = calculate_Vdl_interp(dos_func, Vq_min, Vq_max, C_dl)
        V_dl = V_dl_interp(V_app)
        V_q = V_app - V_dl
        fcn = E->MHC_DOS_integrand(E, λ, V_dl, dos_func; kT=kT, V_q=V_q)
    else
        fcn = E->MHC_DOS_integrand(E, λ, V_app, dos_func; kT=kT, V_q=0)
    end
    quadgk(fcn, E_min, E_max)[1]
end

function MHC_DOS_integrand_redox(E, λ, V_dl, dos_func; ox=true, kT=.026, V_q=0)
    dos_func(E+V_q) * MHC_integrand(E, λ, V_dl, 1 ox; kT=kT)
end




function compute_k_MHC_DOS_redox(E_min, E_max, λ, V_app, dos_func; ox=true, kT=.026, C_q_calc=false, C_dl=10.0, Vq_min=-0.5, Vq_max=0.5)
    if C_q_calc
        V_dl_interp = calculate_Vdl_interp(dos_func, Vq_min, Vq_max, C_dl)
        V_dl = V_dl_interp(V_app)
        V_q = V_app - V_dl
        fcn = E->MHC_DOS_integrand_redox(E, λ, V_dl, dos_func; ox=ox, kT=kT, V_q=V_q)
    else
        fcn = E->MHC_DOS_integrand_redox(E, λ, V_dl, dos_func; ox=ox, kT=kT, V_q=0)
    end
    quadgk(fcn, E_min, E_max)[1]
end

function QC_integrand(E, eVq, dos_func; kT=.026)
    dos_func(E) * sech((E-eVq)/(2*kT))^2/(4*kT)
end

function compute_cq(E_min, E_max, eVq, dos_func; kT=.026)
    fcn = E->QC_integrand(E, eVq, dos_func; kT=kT)
    quadgk(fcn, E_min, E_max)[1]
end


# NOW FOR THE FITTING STUFF

# goal: minimize MSE between experimental and modeled data by choosing best-fit λ
# assumes `exp_data` is an n x 2 array of voltage, current values
# if `dos_func` flag is true, `dos` param will be assumed to be a callable
# otherwise, assumed a float (average dos)
# `p` is a two-element vector with λ plus overall scale
function sq_error(p, E_min, E_max, dos, exp_data, dos_is_func=false)
    # idiot check
    if dos_is_func
        if isempty(methods(dos)) # it's not callable...
            error("It seems you haven't passed a callable DOS. Maybe you meant to set `dos_func`=false?")
        else
            pred_func = V -> p[2]*compute_k_MHC_DOS(E_min, E_max, p[1], V, dos)
        end
    else
        if !isempty(methods(dos[1])) # it's not callable...
            error("It seems you passed a callable DOS when the function expected an average. Maybe you meant to set `dos_func`=true?")
        else
            pred_func = V -> p[2]*compute_k_MHC(E_min, E_max, p[1], V, dos)
        end
    end
    # actually calculate
    V = exp_data[:,1]
    I = exp_data[:,2]
    #sum(log.((pred_func.(V) .- I).^2))
    sum((pred_func.(V) .- I).^2)
end

function best_fit_params(exp_data, dos, E_min, E_max, dos_is_func=false; λ_bounds=(0.01, 0.5), A_bounds=(0.1, 10000))
    V_vals = exp_data[:,1]
    opt_func = p -> sq_error(p, E_min-V_vals[1], E_max-V_vals[end], dos, exp_data, dos_is_func)
    res = bboptimize(opt_func; SearchSpace=RectSearchSpace([λ_bounds, A_bounds]), NumDimensions=2, MaxSteps=100000, TraceInterval=5.0, MinDeltaFitnessTolerance=1e-12)
    best_candidate(res)
end
