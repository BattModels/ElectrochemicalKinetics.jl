using DelimitedFiles
using QuadGK
using Plots
using Statistics
using QuadGK
using Interpolations
include("functions.jl")

for dos_file in string.("DOSes/",readdir("./DOSes"))
    material = dos_file[7:end-7]
    println(material)
    eq_dos, average_eq_dos, min_E_eq, max_E_eq = get_dos(dos_file)

    for Ef in [0.0, 1.0, -1.0, -2.0]
        dos_f, average_dos, min_E, max_E = get_dos(dos_file; Ef=Ef)
        plot_comparison(dos_f, 200, min_E, max_E, plot_title=string("E_f = ", Ef, " eV"))
        savefig(string(material,"_rate_compare_Ef",Ef,".png"))
    end

    # and plot the DOS too
    E_range = range(min_E_eq, max_E_eq-0.01, length=200)
    plot(E_range, eq_dos.dos_interp.(E_range), label="DOS", xlabel="energy", ylabel="DOS")
    savefig(string("figs/", material, "_DOS.png"))
end
