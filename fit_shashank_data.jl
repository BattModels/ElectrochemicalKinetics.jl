using DelimitedFiles
include("functions.jl")

dosfiles = ["DOSes/Li100dos.txt", "DOSes/Li110dos.txt"]
expfiles = string.("shashank_li/", readdir("shashank_li/"))

for expfile in expfiles
    for dosfile in dosfiles
        exp_data = readdlm(expfile, ',')
        dos_f, avg_dos, E_min, E_max = get_dos(dosfile)

        MHC_λ, MHC_A = best_fit_params(exp_data, avg_dos, E_min, E_max, false)
        MHC_DOS_λ, MHC_DOS_A = best_fit_params(exp_data, dos_f, E_min, E_max, true)

        elyte = expfile[13:end-4]
        electrode = dosfile[7:end-7]

        plot_fits(exp_data, dos_f, avg_dos, E_min, E_max, MHC_λ, MHC_A, MHC_DOS_λ, MHC_DOS_A, plot_title=string(elyte, " ", electrode))
        savefig(string("figs/exp_fits/", elyte, "_", electrode, ".png"))
    end
end
