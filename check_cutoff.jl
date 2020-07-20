"""
Do a check to make sure that the bounds of the DOS (provided they're far from Ef) don't affect the final results.
"""

include("functions.jl")

dos_file = "zeeshan_files/relax_pris.dos"
Ef_eq = 0.377 # for now have to pull this manually
E_step = 0.01
kT = .026

cutoff_figprefix = Dict(true => "figs/check_cutoff/cutoffdos/", false => "figs/check_cutoff/fulldos/")

for cut_energy in [true, false]
    eq_dos, average_eq_dos, min_E_eq, max_E_eq = get_dos(dos_file; cut_energy=cut_energy)

    for Ef in [0.0, 1.0, -1.0, -2.0]
        dos_f, average_dos, min_E, max_E = get_dos(dos_file; Ef=Ef, cut_energy=cut_energy)
        println(Ef, " ", min_E, " ", max_E)
        plot_comparison(dos_f, 5, min_E, max_E, average_dos, plot_title=string("E_f = ", Ef, " eV"))
        savefig(string(cutoff_figprefix[cut_energy],"rate_compare_Ef",Ef,".png"))
    end

    # and plot the DOS too
    E_range = range(min_E_eq, max_E_eq-0.01, length=200)
    plot(E_range, eq_dos.(E_range), label="DOS", xlabel="energy", ylabel="DOS")
    savefig(string(cutoff_figprefix[cut_energy], "DOS.png"))
end
