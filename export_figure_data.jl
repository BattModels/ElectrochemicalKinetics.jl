include("analysis_functions.jl")

stepsize=0.005

## FIGURE 1
# Li100 interpolated DOS
dosfile = "Li_100_dos_10smooth.txt"
#dosfile = "DOSes/Li_100_dos.txt"
dos_f, avg_dos, E_min, E_max = get_dos(dosfile)
E_vals = range(E_min, E_max, step=stepsize)
dos_vals = dos_f.(E_vals)
dos_data = hcat(E_vals, dos_vals)
open("figure_data/fig1_Li100/dos_interp.txt", "w") do io
    write(io, "E, dos\n")
    writedlm(io, dos_data)
end;

# MHC fit for ECDEC including Gaussian electrolyte DOSes with fitted parameters
exp_data = readdlm("exp_data/ecdec.txt", ',')
MHC_λ, MHC_A = best_fit_params(exp_data, avg_dos, E_min, E_max, false)
MHC_DOS_λ, MHC_DOS_A = best_fit_params(exp_data, dos_f, E_min, E_max, true)

# save parameter values
params = [MHC_λ, MHC_A, MHC_DOS_λ, MHC_DOS_A]
open("figure_data/fig1_Li100/ecdec_fitparams.txt", "w") do io
    write(io, "MHC_lambda, MHC_A, MHC_DOS_lambda, MHC_DOS_A\n")
    writedlm(io, params)
end;

# calculate and save k vs. V
eη_range = range(-3.5, 3.5, step=stepsize)
MHC_k = [MHC_A * compute_k_MHC(E_min, E_max, MHC_λ, eη, avg_dos) for eη in eη_range]
MHC_DOS_k = [MHC_DOS_A * compute_k_MHC_DOS(E_min, E_max, MHC_DOS_λ, eη, dos_f) for eη in eη_range]

open("figure_data/fig1_Li100/ecdec_fitdata.txt", "w") do io
    write(io, "V, MHC_k, MHC_DOS_k\n")
    writedlm(io, hcat(eη_range, MHC_k, MHC_DOS_k))
end;


## FIGURE 2
# Cu111 DOS
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/Cu111dos.txt")
E_vals = range(E_min, E_max, step=stepsize)
dos_vals = dos_f.(E_vals)
dos_data = hcat(E_vals, dos_vals)
open("figure_data/fig2_Cu111/dos_interp.txt", "w") do io
    write(io, "E, dos\n")
    writedlm(io, dos_data)
end;

# and rate constant comparison, at equilibrium Ef
# for now using same reorg energy from above
V_range = range(-5.0, 5.0, step=stepsize)
MHC_k = [compute_k_MHC(E_min, E_max, MHC_λ, V, avg_dos) for V in V_range]
MHC_DOS_k = [compute_k_MHC_DOS(E_min, E_max, MHC_DOS_λ, V, dos_f) for V in V_range]
open("figure_data/fig2_Cu111/k_vs_V.txt", "w") do io
    write(io, "V, MHC_k, MHC_DOS_k\n")
    writedlm(io, hcat(V_range, MHC_k, MHC_DOS_k))
end;

## FIGURE 3
# SEI (semiconductor) case
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/LiF_Li2CO3dos.txt")
E_vals = range(E_min, E_max, step=stepsize)
dos_vals = dos_f.(E_vals)
dos_data = hcat(E_vals, dos_vals)
open("figure_data/fig3_SEI/dos_interp.txt", "w") do io
    write(io, "E, dos\n")
    writedlm(io, dos_data)
end;

# and rate constant comparison, at equilibrium Ef
# for now using same reorg energy from above
V_range = range(-5.0, 5.0, step=stepsize)
MHC_k = [compute_k_MHC(E_min, E_max, MHC_λ, V, avg_dos) for V in V_range]
MHC_DOS_k = [compute_k_MHC_DOS(E_min, E_max, MHC_DOS_λ, V, dos_f) for V in V_range]
open("figure_data/fig3_SEI/k_vs_V.txt", "w") do io
    write(io, "V, MHC_k, MHC_DOS_k\n")
    writedlm(io, hcat(V_range, MHC_k, MHC_DOS_k))
end;

## FIGURE 4
# impact of changing eq. Ef on SEI
# first do DOS again
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/LiF_Li2CO3dos.txt"; Ef=-2.0)
E_vals = range(E_min, E_max, step=stepsize)
dos_vals = dos_f.(E_vals)
dos_data = hcat(E_vals, dos_vals)
open("figure_data/fig4_params/SEI_dos_interp_Ef-2.txt", "w") do io
    write(io, "E, dos\n")
    writedlm(io, dos_data)
end;

# then recompute k's
MHC_k = [compute_k_MHC(E_min, E_max, MHC_λ, V, avg_dos) for V in V_range]
MHC_DOS_k = [compute_k_MHC_DOS(E_min, E_max, MHC_DOS_λ, V, dos_f) for V in V_range]
open("figure_data/fig4_params/SEI_k_vs_V_Ef-2.txt", "w") do io
    write(io, "V, MHC_k, MHC_DOS_k\n")
    writedlm(io, hcat(V_range, MHC_k, MHC_DOS_k))
end;


# and now, impact of changing λ on the Cu111 case
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/Cu111dos.txt")
V_range = range(-5.0, 5.0, step=stepsize)
MHC_k = [compute_k_MHC(E_min, E_max, MHC_DOS_λ*2, V, avg_dos) for V in V_range]
MHC_DOS_k = [compute_k_MHC_DOS(E_min, E_max, MHC_DOS_λ*2, V, dos_f) for V in V_range]
open("figure_data/fig5/Cu111_k_vs_V_doublelambda.txt", "w") do io
    write(io, "V, MHC_k, MHC_DOS_k\n")
    writedlm(io, hcat(V_range, MHC_k, MHC_DOS_k))
end;
