include("functions.jl")

for dos_file in string.("DOSes/",readdir("./DOSes"))
    material = dos_file[7:end-7]
    println(material)
    if occursin("Li3Ag", dos_file)
        eq_dos, average_eq_dos, min_E_eq, max_E_eq = get_dos(dos_file, cut_energy=true)
    else
        eq_dos, average_eq_dos, min_E_eq, max_E_eq = get_dos(dos_file)
    end

    for Ef in [0.0, 1.0, -1.0, -2.0]
        dos_f, average_dos, min_E, max_E = get_dos(dos_file; Ef=Ef)
        plot_comparison(dos_f, 5, min_E, max_E, average_dos, plot_title=string("E_f = ", Ef, " eV"))
        savefig(string("figs/rate_compare/", material,"_Ef",Ef,".png"))
    end

    # and plot the DOS too
    E_range = range(min_E_eq, max_E_eq-0.01, length=200)
    plot(E_range, eq_dos.(E_range), label="DOS", xlabel="energy", ylabel="DOS")
    savefig(string("figs/DOS/", material, ".png"))
end
