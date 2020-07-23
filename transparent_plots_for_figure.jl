using LaTeXStrings
include("functions.jl")

# make transparent plots with no axes,
# just the lines, to pull over for a schematic figure
# keywords: background_color=:transparent, axis=nothing, border=:none, legend=false

min_E = -5
max_E = 5
E_range = range(min_E, max_E, length=500)
kT = 0.1

plot_args = Dict(:background_color=>:transparent, :axis=>nothing, :border=>:none, :legend=>false, :size=>(500,100))

# Fermi-Dirac
fd_vals_e = fermi_dirac.(E_range; kT=kT)
fd_vals_h = 1 .- fd_vals_e
plot(E_range, fd_vals_e; plot_args...)
savefig("figs/for_schematic/fd.png")
plot(E_range, hcat(fd_vals_e, fd_vals_h); plot_args...)
savefig("figs/for_schematic/fd_both.png")

# DOS for semiconductor
dos_file = "zeeshan_files/relax_pris.dos"
dos_f_sc, average_dos, min_E, max_E = get_dos(dos_file; Ef=0.377)
dos_vals = dos_f_sc.(E_range)
plot(E_range, dos_vals; color=:black, plot_args...)
savefig("figs/for_schematic/semi_dos.png")
plot(E_range, hcat(dos_vals.*fd_vals_e, dos_vals.*fd_vals_h); plot_args...)
savefig("figs/for_schematic/semi_dos_fdmult.png")

# DOS for metal
dos_file = "DOSes/Cu100dos.txt"
dos_f_m, average_dos, min_E, max_E = get_dos(dos_file)
dos_vals = dos_f_m.(E_range)
plot(E_range, dos_vals; color=:green, plot_args...)
savefig("figs/for_schematic/metal_dos.png")
plot(E_range, hcat(dos_vals.*fd_vals_e, dos_vals.*fd_vals_h); plot_args...)
savefig("figs/for_schematic/metal_dos_fdmult.png")

# Gaussian electrolyte DOSes
λ = 0.2
ox_dos_vals_eq = marcus_integrand.(E_range, λ, 0, true; kT=kT)
red_dos_vals_eq = marcus_integrand.(E_range, λ, 0, false; kT=kT)
ox_dos_vals = marcus_integrand.(E_range, λ, -2.0, true; kT=kT)
red_dos_vals = marcus_integrand.(E_range, λ, -2.0, false; kT=kT)
plot(E_range, hcat(red_dos_vals_eq, ox_dos_vals_eq); plot_args...)
savefig("figs/for_schematic/elyte_dos_eq.png")
plot(E_range, hcat(red_dos_vals, ox_dos_vals); linestyle=:dot, plot_args...)
savefig("figs/for_schematic/elyte_dos_wpot.png")

# okay now plot k vs. eη at the two different voltages...
# first for semiconductor...
dos_file = "zeeshan_files/relax_pris.dos"
dos_f_sc, average_dos, min_E, max_E = get_dos(dos_file; Ef=0.377)
dos_f_sc_noneq, average_dos, min_E, max_E = get_dos(dos_file; Ef=0.377-2.0)

min_eη = -5
max_eη = 5
eη_range = range(min_eη, max_eη, length=500)
k_eq = [compute_k_MHC_DOS(min_eη, max_eη, λ, eη, dos_f_sc) for eη in eη_range]
k_noneq = [compute_k_MHC_DOS(min_eη, max_eη, λ, eη, dos_f_sc_noneq) for eη in eη_range]

plot(eη_range, hcat(log10.(abs.(k_eq)), log10.(abs.(k_noneq))), label=[L"E_f(eq)" L"E_f(eq)-2.0V"], color=[:green :purple], legend=:bottomright, xlabel="eη", ylabel="log(k)")

# then for metal...
dos_file = "DOSes/Cu100dos.txt"
dos_f_m, average_dos, min_E, max_E = get_dos(dos_file)
dos_f_m_noneq, avg, min, max = get_dos(dos_file; Ef=-2.0)
min_eη = -5
max_eη = 5
eη_range = range(min_eη, max_eη, length=500)
k_eq = [compute_k_MHC_DOS(min_eη, max_eη, λ, eη, dos_f_m) for eη in eη_range]
k_noneq = [compute_k_MHC_DOS(min_eη, max_eη, λ, eη, dos_f_m_noneq) for eη in eη_range]
plot(eη_range, hcat(log10.(abs.(k_eq)), log10.(abs.(k_noneq))), label=[L"E_f(eq)" L"E_f(eq)-2.0V"], color=[:green :purple], legend=:bottomright, xlabel="eη", ylabel="log(k)")
