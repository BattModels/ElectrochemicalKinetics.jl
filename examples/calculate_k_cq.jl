include("analysis_fcns.jl")

stepsize = 0.001
lambda = 0.82
V_app = 0.4
r = 0.45

## bernal stacked graphene
dos_f, avg_dos, E_min, E_max = get_dos("../data/DOSes/bernal_graphene.txt"; Ef=0)
k = compute_k_MHC_DOS(-r, r, lambda, V_app, dos_f; kT=.026,  C_q_calc=true)

println("Rate constant: ", k)
