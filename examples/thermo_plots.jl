# basic setup

include("thermo.jl")

xs = 0.01:0.01:0.95

bv =  ButlerVolmer(20)
m = Marcus(350, 0.3)
amhc = AsymptoticMarcusHushChidsey(5000, 0.3)

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


# plot multiple chemical potentials at the same current
I = 20

μ_bv = μ_kinetic(I, bv)
μ_m = μ_kinetic(I, m)
μ_amhc = μ_kinetic(I, amhc)

plot(xs, [μ_bv(xs) μ_m(xs) μ_amhc(xs)], label=["Butler-Volmer" "Marcus" "Marcus-Hush-Chidsey"], legend=:topleft)
xlabel!("x")
ylabel!("μ")

# multiple currents for one model type
I_vals = [10, 20, 30, 40]
μ_fcns = Dict(I => μ_kinetic(I, bv) for I in I_vals)
plot(xs, [μ_fcns[I](xs) for I in I_vals], label=I_vals', legend=:topleft, legendtitle="I")
xlabel!("x")
ylabel!("μ")

# make phase diagrams!  