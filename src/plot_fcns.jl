using Plots
include("analysis_functions.jl")

# plotting shortcut to make MHC and MHC+DOS curves
function plot_comparison(dos_func, eη_max, E_min, E_max, average_dos; kT=.026, λ=.26, length=500, plot_title="")
    eη_range = range(-eη_max, eη_max, length=length)
    MHC_k = [compute_k_MHC(E_min, E_max, λ, eη, average_dos; kT=kT) for eη in eη_range]
    MHC_DOS_k = [compute_k_MHC_DOS(E_min, E_max, λ, eη, dos_func; kT=kT) for eη in eη_range]
    if any(MHC_k.<0) | any(MHC_DOS_k.<0)
        println("negative k...uh-oh...")
        println(sum(MHC_k.<0))
        println(sum(MHC_DOS_k.<0))
    end
    plot(eη_range, hcat(log10.(abs.(MHC_k)), log10.(abs.(MHC_DOS_k))), label=["MHC" "MHC+DOS"], legend=:bottomright, xlabel="eη", ylabel="log(k)", title=plot_title, color=[:purple :green])
end

# if you want to visualize each term in integrand separately
function plot_integrand_components(dos_func, E_min, E_max, average_dos; eη=0, kT=.026, λ=.26, length=500, plot_title="")
    E_range = range(E_min, E_max, length=length)
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
    plot(xs, ys, seriestype=[:scatter :line :line ], label=["experiment" string("MHC: λ=",round(MHC_λ, digits=3), "; A=", round(MHC_A, digits=2)) string("MHC+DOS: λ=", round(MHC_DOS_λ, digits=3), "; MHC_DOS_A=", round(MHC_DOS_A, digits=2))], xlabel="V", ylabel="log(k or I)", yscale=:log10, leg=:bottomright, title=plot_title)
end
