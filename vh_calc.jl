include("analysis_fcns.jl")

stepsize=0.001

Vq_range = range(-0.3, 0.3, step=0.001)


# bernal (AB) graphene bilayer
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/bernal_graphene.txt"; Ef=0)
CQ = [compute_cq(E_min, E_max, V, dos_f)*1.602*10 for V in Vq_range]
Vh = []
Vappl = []
for i = 1:length(Vq_range)
    push!(Vappl, Vq_range[i]*((CQ[i]/10.0) + 1))
end
for i = 1:length(Vq_range)
    push!(Vh, Vappl[i] - Vq_range[i])
end
open("figure_data/bernal_graphene_vh_vs_vappl.txt", "w") do io
    write(io, "Vappl, Vh\n")
    writedlm(io, hcat(Vappl, Vh))
end;
