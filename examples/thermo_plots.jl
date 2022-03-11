# basic setup
using Plots.PlotMeasures
include("thermo.jl")

xs = 0.01:0.01:0.95

bv =  ButlerVolmer(20)
m = Marcus(350, 0.3)
amhc = AsymptoticMarcusHushChidsey(5000, 0.3)

## CURRENT VS. OVERPOTENTIAL

# plot wide range of overpotentials
Vs = -0.5:0.0051:0.5
plot(Vs, [bv(Vs), m(Vs), amhc(Vs)], yaxis=:log, label=["Butler-Volmer" "Marcus" "Marcus-Hush-Chidsey"], size=(400, 400), lw=3, xticks=[], yticks=[])
ylims!((1, 1e6))
xlabel!("η [V]")
ylabel!("log(I)")

# narrower range of overpotentials
Vs = 0.005:0.005:0.3
plot(Vs, [bv(Vs), m(Vs), amhc(Vs)], yaxis=:log, label=["Butler-Volmer" "Marcus" "Marcus-Hush-Chidsey"], legend=:topleft, size=(400, 400), lw=3, xticks=[], yticks=[])
xlabel!("η [V]")
ylabel!("log(I)")


# THERMODYNAMIC FUNCTIONS
# TODO: fix colors so reference case is always blue

I = 20

μ_bv = μ_kinetic(I, bv)
μ_m = μ_kinetic(I, m)
μ_amhc = μ_kinetic(I, amhc)

μ_vals_bv = μ_bv(xs)
μ_vals_m = μ_m(xs)
μ_vals_amhc = μ_amhc(xs)

ref_color = palette(:tab10)[1]

# compare thermodynamic and kinetic cases
plot(xs, [μ_thermo(xs), μ_vals_bv .- μ_thermo(xs), μ_vals_bv], 
    lw=3, 
    legend=:topleft,
    xlabel="x", 
    label = ["Thermodynamic μ" "Δϕ" "Kinetic μ (Butler-Volmer)"], 
    linecolor = [:grey :darkorchid4 ref_color],
    yticks=[], ylims=[-0.02, 0.3],
)

# plot multiple chemical potentials at the same current
plot(xs, [μ_vals_m, μ_vals_amhc, μ_vals_bv], 
    label=["Marcus" "Marcus-Hush-Chidsey" "Butler-Volmer"], legend=:topleft, 
    lw=3,  
    xlabel="x", ylabel="μ", 
    yticks=[], ylims=[-0.02,0.3],
    linecolor = [palette(:tab10)[2] palette(:tab10)[3] ref_color],
    legendtitle="Vary model at I=I₀", 
    margin=3mm,
)

# multiple currents for one model type
I_vals = [20, 100, 10]
μ_fcns = Dict(I => μ_kinetic(I, bv) for I in I_vals)
plot(xs, [μ_fcns[I](xs) for I in I_vals], 
    label=["I₀" "10I₀" "0.5I₀"],
    legend=:topleft, legendtitle="Butler-Volmer, vary I", xlabel="x", ylabel="μ", 
    yticks=[], 
    linecolor = [ref_color :maroon :goldenrod2],
    lw=3, 
    ylims=[-0.02,0.3],
    margin=3mm
)

## PHASE DIAGRAM - BUTLER-VOLMER
# first a 2D BV one

I_vals = 2.5:2.5:75
# I_vals = [2.5, 20, 75]
pb1 = []
pb2 = []
for I in I_vals
    println(I)
    pbs = find_phase_boundaries(I, bv, T=330)
    push!(pb1, pbs[1])
    push!(pb2, pbs[2])
end

open("pbs_bv_T330.txt", "w") do io
    writedlm(io, [I_vals pb1 pb2])
end


I_vals = 2.5:2.5:100
pb1 = []
pb2 = []
for I in I_vals
    println(I)
    pbs = find_phase_boundaries(I, bv, T=310)
    push!(pb1, pbs[1])
    push!(pb2, pbs[2])
end

open("pbs_bv_T310.txt", "w") do io
    writedlm(io, [I_vals pb1 pb2])
end

# I_vals = 2.5:2.5:100
# pb1 = []
# pb2 = []
# for I in I_vals
#     println(I)
#     pbs = find_phase_boundaries(I, bv, T=320)
#     push!(pb1, pbs[1])
#     push!(pb2, pbs[2])
# end

# open("pbs_bv_T320.txt", "w") do io
#     writedlm(io, [I_vals pb1 pb2])
# end

# I_vals = 2:2:30
# for T in 350:10:430
#     pb1 = []
#     pb2 = []
#     for I in I_vals
#         println(T, I)
#         pbs = find_phase_boundaries(I, bv, T=T)
#         push!(pb1, pbs[1])
#         push!(pb2, pbs[2])
#     end

#     open("pbs_bv_T$T.txt", "w") do io
#         writedlm(io, [I_vals pb1 pb2])
#     end
# end

# all_data = []
# for T in 310:10:400    
#     open("pbs_bv_T$T.txt", "r") do io
#         data = readdlm(io)
#         data = hcat(T .* ones(size(data,1)), data)
#         push!(all_data, data)
#     end
# end

# all_data = vcat(all_data...)
# xy = all_data[:, 1:2]
# z1 = all_data[:, 3]
# z2 = all_data[:, 4]

# all_data = vcat([xy z1], [xy z2])

open("pbs_bv.txt", "w") do io
    writedlm(io, all_data)
end

data = readdlm("pbs_bv.txt")

# mesh3d(all_data[:,1], all_data[:,3], all_data[:,2], xlabel="T", ylabel="x", zlabel="I", zticks=[], legend=:none)
Plots.plot(data[:,1], data[:,3], data[:,2], st = :scatter3d, xlabel="T", ylabel="x", zlabel="I", legend=:none)

## PHASE DIAGRAM - MARCUS
I_vals = 4:1:16
pb1 = []
pb2 = []
for I in I_vals
    println(I)
    pbs = find_phase_boundaries(I, m, T=370)
    push!(pb1, pbs[1])
    push!(pb2, pbs[2])
end

open("pbs_m_T370.txt", "w") do io
    writedlm(io, [I_vals pb1 pb2])
end