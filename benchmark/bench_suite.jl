using BenchmarkTools
using ElectrochemicalKinetics

const suite = BenchmarkGroup()
time_mult = 1

suite["integrals"] = BenchmarkGroup()
suite["integrals"]["scale"] = @benchmarkable ElectrochemicalKinetics.scale(0.1, 0.2) seconds=time_mult

# reusing the same ones from phase_diagram_tests
bv = ButlerVolmer(300)
m = Marcus(5000, 0.3)
amhc = AsymptoticMarcusHushChidsey(70000, 0.3)
ni_kms = [bv, m, amhc]

# rate constant benchmarks
# TODO: add vector of models cases
suite["rate constants"] = BenchmarkGroup()
suite["rate constants"]["non-integral models"] = BenchmarkGroup()
suite["rate constants"]["non-integral models"]["scalar V"] = BenchmarkGroup()
suite["rate constants"]["non-integral models"]["vector V"] = BenchmarkGroup()

for km in ni_kms
    modeltype = typeof(km)
    suite["rate constants"]["non-integral models"]["scalar V"][modeltype] = @benchmarkable compute_k(0.1, $km) seconds=time_mult*0.1
    suite["rate constants"]["non-integral models"]["vector V"][modeltype] = @benchmarkable compute_k(0.01:0.01:0.2, $km) seconds=time_mult
end

mhc = MarcusHushChidsey(70000, 0.3)
mhcd_flatdos = MarcusHushChidseyDOS(70000, 0.3, [-5 1; 5 1])
mhcd_Cu = MarcusHushChidseyDOS(10000, 0.3, string(@__DIR__, "/../data/DOSes/Cu_111_dos.txt"))
i_kms = [mhc, mhcd_flatdos, mhcd_Cu]
names = Dict(mhc=>"MarcusHushChidsey", mhcd_flatdos=>"MarcusHushChidseyDOS: flat dos", mhcd_Cu=>"MarcusHushChidseyDOS: Cu111")

suite["rate constants"]["integral models"] = BenchmarkGroup()
suite["rate constants"]["integral models"]["scalar V"] = BenchmarkGroup()
suite["rate constants"]["integral models"]["vector V"] = BenchmarkGroup()

for km in i_kms
    modeltype = names[km]
    suite["rate constants"]["integral models"]["scalar V"][modeltype] = @benchmarkable compute_k(0.1, $km) seconds=time_mult*8
    suite["rate constants"]["integral models"]["vector V"][modeltype] = @benchmarkable compute_k(0.01:0.01:0.2, $km) seconds=time_mult*11
end

# fit_overpotential benchmarks
# TODO: add vector cases
k = 2e3

suite["overpotential fitting"] = BenchmarkGroup()
suite["overpotential fitting"]["non-integral models"] = BenchmarkGroup()
suite["overpotential fitting"]["non-integral models"]["non-AD"] = BenchmarkGroup()
suite["overpotential fitting"]["non-integral models"]["AD"] = BenchmarkGroup()

for km in ni_kms
    modeltype = typeof(km)
    suite["overpotential fitting"]["non-integral models"]["non-AD"][modeltype] = @benchmarkable fit_overpotential($km, $k, autodiff=false) seconds=time_mult*3
    suite["overpotential fitting"]["non-integral models"]["AD"][modeltype] = @benchmarkable fit_overpotential($km, $k, autodiff=true) seconds=time_mult*5
end

suite["overpotential fitting"]["integral models"] = BenchmarkGroup()
suite["overpotential fitting"]["integral models"]["non-AD"] = BenchmarkGroup()
suite["overpotential fitting"]["integral models"]["AD"] = BenchmarkGroup()
for km in i_kms[2:end] # skipping MHC for now
    modeltype = names[km]
    suite["overpotential fitting"]["integral models"]["non-AD"][modeltype] = @benchmarkable fit_overpotential($km, $k, autodiff=false) seconds=time_mult*10
    suite["overpotential fitting"]["integral models"]["AD"][modeltype] = @benchmarkable fit_overpotential($km, $k, autodiff=true) seconds=time_mult*20
end
