using ElectrochemicalKinetics

stepsize = 0.001
λ = 0.82
V_app = 0.4
r = 0.45

## bernal stacked graphene
dos_file = "data/DOSes/bernal_graphene.txt"
model = MarcusHushChidseyDOS(λ, dos_file)
k = compute_k_cq(-r, r, V_app, model)

println("Rate constant: ", k)
