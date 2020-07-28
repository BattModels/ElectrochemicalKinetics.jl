using Interact
using Plots
using LaTeXStrings
include("functions.jl")

dosfiles = string.("DOSes/", readdir("./DOSes/"))
electrode_surfaces = [dosfile[7:end-7] for dosfile in dosfiles]

@manipulate for Ef_eq=-2.0:0.05:1.0, λ=0.05:0.01:0.5, T=100:20:400, electrode_surface=electrode_surfaces, show_electrolyte_DOS=true
    kT = 8.617e-5*T
    dosfile = string("DOSes/", electrode_surface, "dos.txt")
    dos_f, average_dos, min_E, max_E = get_dos(dosfile; Ef=Ef_eq)
    if occursin("Li1",electrode_surface)
        E_range=1.5
    else
        E_range=4
    end
    E_vals = range(-E_range, E_range, length=500)
    dos_vals=dos_f.(E_vals)
    fd_vals_e = fermi_dirac.(E_vals; kT=kT)
    fd_vals_h = 1 .- fd_vals_e
    dos_vals_e = fd_vals_e .* dos_vals
    dos_vals_h = fd_vals_h .* dos_vals
    if show_electrolyte_DOS
        ox_dos_vals = 0.5.*maximum(dos_vals).*marcus_integrand.(E_vals, λ, 0, true; kT=kT)
        red_dos_vals = 0.5.*maximum(dos_vals).*marcus_integrand.(E_vals, λ, 0, false; kT=kT)
        p1 = plot(E_vals, hcat(dos_vals_e, dos_vals_h, ox_dos_vals, red_dos_vals), label=["electrons" "holes" "reduction" "oxidation"], title="Electrode state occupations", xlabel="E-E_f", ylabel="# (arb.)", color=[:blue :orange :blue :orange], linestyle=[:solid :solid :dash :dash])
    else
        p1 = plot(E_vals, hcat(dos_vals_e, dos_vals_h), label=["electrons" "holes"], title="Electrode state occupations", xlabel=L"E-E_f", ylabel="# (arb.)")
    end
    p3 = plot_comparison(dos_f, E_range, min_E, max_E, average_dos; kT=kT, λ=λ, plot_title="Rate constants")
    plot(p1, p3, layout=(2,1), size=(500,500))
end

#w = Window()
#body!(w, ui);
