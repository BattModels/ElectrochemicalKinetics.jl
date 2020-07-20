using DelimitedFiles
using QuadGK
using Plots
using Statistics
using QuadGK
using Interpolations
include("functions.jl")

kT = 0.026
λ = 10

# semiconductor case, equilibrium
dos_file = "zeeshan_files/relax_pris.dos"
eq_dos, average_eq_dos, min_E_eq, max_E_eq = get_dos(dos_file; Ef=0.377)

E_range = range(-5, 5, length=200)
x_range = E_range ./kT
dos_vals = eq_dos.dos_interp.(E_range)

fd_vals = fermi_dirac.(x_range)
ox_dos_vals = marcus_integrand.(x_range, λ, 0, true)
red_dos_vals = marcus_integrand.(x_range, λ, 0, false)

to_plot = hcat(dos_vals, fd_vals, ox_dos_vals, red_dos_vals)
plot(x_range, to_plot, layout=(4,1))
