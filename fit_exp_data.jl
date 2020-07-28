using DelimitedFiles
using Statistics
include("functions.jl")

dosfiles = ["DOSes/Li_100_dos.txt"]
expfiles = string.("exp_data/", readdir("exp_data/"))

for expfile in expfiles
    for dosfile in dosfiles
        exp_data = readdlm(expfile, ',')
        dos_f, avg_dos, E_min, E_max = get_dos(dosfile)

        MHC_λ, MHC_A = best_fit_params(exp_data, avg_dos, E_min, E_max, false)
        MHC_DOS_λ, MHC_DOS_A = best_fit_params(exp_data, dos_f, E_min, E_max, true)

        # calculate RMSE
        V_vals = exp_data[:,1]
        I = exp_data[:,2]
        MHC_pred = [MHC_A*compute_k_MHC(E_min, E_max, MHC_λ, V, avg_dos) for V in V_vals]
        MHC_DOS_pred = [MHC_DOS_A*compute_k_MHC_DOS(E_min, E_max, MHC_DOS_λ, V, dos_f) for V in V_vals]
        MHC_rmse = sqrt(mean((MHC_pred .- I).^2))
        MHC_DOS_rmse = sqrt(mean((MHC_DOS_pred .- I).^2))

        elyte = expfile[13:end-4]
        electrode = dosfile[7:end-7]
        plot_fits(exp_data, dos_f, avg_dos, E_min, E_max, MHC_λ, MHC_A, MHC_DOS_λ, MHC_DOS_A, plot_title=string(elyte, " ", electrode, "; MHC RMSE=", round(MHC_rmse, digits=2), "; MHC_DOS RMSE=", round(MHC_DOS_rmse, digits=2)))
        savefig(string("figs/exp_fits/", elyte, "_", electrode, ".png"))
    end
end
