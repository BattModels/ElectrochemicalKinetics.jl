using Statistics
using Interpolations
using QuadGK
using Plots
using LaTeXStrings
using DelimitedFiles
using BlackBoxOptim

#TODO: clean up max/min variable names

# assumes `dos_file` has two columns: energy (relative to equilibrium Ef) and DOS
function get_dos(dos_file; Ef=0, cut_energy=false)
    dos_data = readdlm(dos_file, Float64, skipstart=1)
    # recenter so Ef=0
    dos_data[:,1] = dos_data[:,1] .- Ef
    max_E = dos_data[end,1]
    if cut_energy
        # for now, we'll pick a symmetric energy range about 0 (i.e. the new Ef)
        len_keep = sum(dos_data[:,1] .> -max_E)
        # cut off
        dos_data = dos_data[end-len_keep:end,:]
    end
    min_E = dos_data[1,1]
    average_dos = mean(dos_data[:,2]) # for whole structure
    # interpolate
    E_step = mean(dos_data[2:end,1].-dos_data[1:end-1,1])
    dos_interp = scale(interpolate(dos_data[:,2], BSpline(Linear())), range(min_E, max_E+0.0001, step=E_step))
    return dos_interp, average_dos, min_E, max_E
end

# assuming x is in units of kT
function fermi_dirac(E; kT=.026)
    1 / (1 + exp(E/kT))
end

# all arguments in units of kT
function marcus_integrand(E, λ, eη, ox=true; kT=.026)
    if ox # oxidative direction
        arg = -(λ-eη+E)^2 / (4*λ*kT)
    else # reductive direction
        arg = -(λ+eη-E)^2 / (4*λ*kT)
    end
    exp(arg)
end

function MHC_integrand(E, λ, eη, average_dos; kT=.026)
    marcus_ox = marcus_integrand(E, λ, eη, true; kT=kT)
    marcus_red = marcus_integrand(E, λ, eη, false; kT=kT)
    fd_ox = 1-fermi_dirac(E; kT=kT) # holes
    fd_red = fermi_dirac(E; kT=kT) # electrons
    sign(eη) * average_dos * (marcus_ox*fd_ox - marcus_red*fd_red)
end

function compute_k_MHC(E_min, E_max, λ, eη, average_dos; kT=.026)
    fcn = E->MHC_integrand(E, λ, eη, average_dos; kT=.026)
    quadgk(fcn, E_min, E_max)[1] # return format is (value, error_bound)
end

# compute with DOS by evaluating the function we already have for a unit DOS and multiplying
# by the interpolated value
function MHC_DOS_integrand(E, λ, eη, dos_func; kT=.026)
    dos_func(E) * MHC_integrand(E, λ, eη, 1; kT=kT)
end

function compute_k_MHC_DOS(E_min, E_max, λ, eη, dos_func; kT=.026)
    fcn = E->MHC_DOS_integrand(E, λ, eη, dos_func; kT=kT)
    quadgk(fcn, E_min, E_max)[1]
end

function plot_comparison(dos_func, eη_max, min_E, max_E, average_dos; kT=.026, λ=.26, length=500, plot_title="")
    eη_range = range(-eη_max, eη_max, length=length)
    MHC_k = [compute_k_MHC(min_E, max_E, λ, eη, average_dos; kT=kT) for eη in eη_range]
    MHC_DOS_k = [compute_k_MHC_DOS(min_E, max_E, λ, eη, dos_func; kT=kT) for eη in eη_range]
    if any(MHC_k.<0) | any(MHC_DOS_k.<0)
        println("negative k...uh-oh...")
        println(sum(MHC_k.<0))
        println(sum(MHC_DOS_k.<0))
    end
    plot(eη_range, hcat(log10.(abs.(MHC_k)), log10.(abs.(MHC_DOS_k))), label=["MHC" "MHC+DOS"], legend=:bottomright, xlabel="eη", ylabel="log(k)", title=plot_title)
end

function plot_integrand_components(dos_func, min_E, max_E, average_dos; eη=0, kT=.026, λ=.26, length=500, plot_title="")
    E_range = range(min_E, max_E, length=length)
    dos_vals = dos_func.(E_range)
    p1 = plot(E_range, dos_vals, color=:black, ylabel="electrode states", label="electrons/holes", title=plot_title)

    fd_vals_e = fermi_dirac.(E_range; kT=kT)
    fd_vals_h = 1 .- fd_vals_e
    p2 = plot(E_range, hcat(fd_vals_e, fd_vals_h), label=["electrons" "holes"], ylabel="occupation", legend=:right)

    ox_dos_vals = marcus_integrand.(E_range, λ, eη, true; kT=kT)
    red_dos_vals = marcus_integrand.(E_range, λ, eη, false; kT=kT)
    p3 = plot(E_range, hcat(red_dos_vals, ox_dos_vals), label=["reduction" "oxidation"], ylabel="electrolyte states")

    plot(p1, p2, p3, layout=(3,1))
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

function plot_fits(exp_data, dos_f, avg_dos, E_min, E_max, MHC_λ, MHC_A, MHC_DOS_λ, MHC_DOS_A; plot_title="")
    # compute k's
    V = exp_data[:,1]
    V_mag = 1.1 * maximum(abs.(V))
    V_range = range(-V_mag, V_mag, length=200)
    MHC_k = [MHC_A*compute_k_MHC(E_min, E_max, MHC_λ, V, avg_dos) for V in V_range]
    MHC_DOS_k = [MHC_DOS_A*compute_k_MHC_DOS(E_min, E_max, MHC_DOS_λ, V, dos_f) for V in V_range]

    # scatter plot of experimental data, lines for fits
    xs = Vector[V, V_range, V_range]
    ys = Vector[exp_data[:,2], MHC_k, MHC_DOS_k]
    #plot(V, log10.(exp_data[:,2]), seriestype=:scatter, label="experiment")
    #plot!(V_range, hcat(log10.(abs.(MHC_k)), log10.(abs.(MHC_DOS_k))), label=[string(L"MHC: $\lambda=$",MHC_λ) string(L"MHC+DOS: $\lambda=$", MHC_DOS_λ)], xlabel="V", ylabel="log(k or I)")
    plot(xs, ys, seriestype=[:scatter :line :line ], label=["experiment" string("MHC: λ=",round(MHC_λ, digits=3)) string("MHC+DOS: λ=", round(MHC_DOS_λ, digits=3))], xlabel="V", ylabel="log(k or I)", yscale=:log10, leg=:bottomright, title=plot_title)
end
