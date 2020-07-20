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
    # function to get the interpolated value for the nondimensional argument (need to assume a temperature)
    dos(x, kT=0.026) = dos_interp(kT*x)
    return dos, average_dos, min_E, max_E
end

# assuming x is in units of kT
function fermi_dirac(x)
    1 / (1 + exp(x))
end

# all arguments in units of kT
function marcus_integrand(x, λ, eη, ox=true)
    if ox # oxidative direction
        arg = -(λ-eη+x)^2 / (4*λ)
    else # reductive direction
        arg = -(λ+eη-x)^2 / (4*λ)
    end
    exp(arg)
end

function MHC_integrand(x, λ, eη, average_dos)
    marcus_ox = marcus_integrand(x, λ, eη, true)
    marcus_red = marcus_integrand(x, λ, eη, false)
    fd_ox = 1-fermi_dirac(x) # holes
    fd_red = fermi_dirac(x) # electrons
    sign(eη) * average_dos * (marcus_ox*fd_ox - marcus_red*fd_red)
end

function compute_k_MHC(x_min, x_max, λ, eη, average_dos)
    fcn = x->MHC_integrand(x, λ, eη, average_dos)
    quadgk(fcn, x_min, x_max)[1] # return format is (value, error_bound)
end

# compute with DOS by evaluating the function we already have for a unit DOS and multiplying
# by the interpolated value
function MHC_DOS_integrand(dos_func, x, λ, eη, kT=.026)
    dos_func(x, kT) * MHC_integrand(x, λ, eη, 1)
end

function compute_k_MHC_DOS(dos_func, x_min, x_max, λ, eη, kT=.026)
    fcn = x->MHC_DOS_integrand(dos_func, x, λ, eη, kT)
    quadgk(fcn, x_min, x_max)[1]
end

function plot_comparison(dos_func, eη_max, min_E, max_E, average_dos; kT=.026, length=500, plot_title="")
    eη_range = range(-eη_max, eη_max, length=length)
    MHC_k = [compute_k_MHC(min_E/kT, max_E/kT, 10, eη, average_dos) for eη in eη_range]
    MHC_DOS_k = [compute_k_MHC_DOS(dos_func, min_E/kT, max_E/kT, 10, eη) for eη in eη_range]
    if any(MHC_k.<0) | any(MHC_DOS_k.<0)
        println("negative k...uh-oh...")
        println(sum(MHC_k.<0))
        println(sum(MHC_DOS_k.<0))
    end
    plot(eη_range, hcat(log10.(abs.(MHC_k)), log10.(abs.(MHC_DOS_k))), label=["MHC" "MHC+DOS"], legend=:bottomright, xlabel="eη", ylabel="log(k)", title=plot_title)
end
