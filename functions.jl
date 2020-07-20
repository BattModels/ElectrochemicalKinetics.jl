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
function MHC_DOS_integrand(dos_func, E, λ, eη; kT=.026)
    dos_func(E) * MHC_integrand(E, λ, eη, 1; kT=kT)
end

function compute_k_MHC_DOS(dos_func, E_min, E_max, λ, eη; kT=.026)
    fcn = E->MHC_DOS_integrand(dos_func, E, λ, eη; kT=kT)
    quadgk(fcn, E_min, E_max)[1]
end

function plot_comparison(dos_func, eη_max, min_E, max_E, average_dos; kT=.026, λ=.26, length=500, plot_title="")
    eη_range = range(-eη_max, eη_max, length=length)
    MHC_k = [compute_k_MHC(min_E, max_E, λ, eη, average_dos) for eη in eη_range]
    MHC_DOS_k = [compute_k_MHC_DOS(dos_func, min_E, max_E, λ, eη) for eη in eη_range]
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
