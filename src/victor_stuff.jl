using DelimitedFiles
using ElectrochemicalKinetics
using Plots

dos_data = readdlm("/Users/rkurchin/Downloads/dos.txt")
E = dos_data[:,1]

dos = DOSData("/Users/rkurchin/Downloads/dos.txt")
λ = .5
A = 3.4376761519
model = MarcusHushChidseyDOS(A, λ, dos)

V = -3.0:0.05:3.0
k_vals = compute_k(V, model, kT=kT)
plot(V, k_vals)

η_vals = [-0.2, 0.0, 0.2]

kT = .0258

marcus_term(E) = exp.(-(( (mhcd.λ .- V_dl) .+ E) .^ 2) ./ (4 * mhcd.λ * kT))
marcus_factors = marcus_term.(E)

dos_term = (1 .- fermi_dirac(E; kT=kT)) .* model.dos.interp_func.(E)
