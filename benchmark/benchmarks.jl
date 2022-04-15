using BenchmarkTools
using ElectrochemicalKinetics

const SUITE = BenchmarkGroup()

SUITE["integrals"] = BenchmarkGroup()
SUITE["integrals"]["scale"] = @benchmarkable ElectrochemicalKinetics.scale(0.1, 0.2)

# reusing the same ones from phase_diagram_tests
bv = ButlerVolmer(300)
m = Marcus(5000, 0.3)
amhc = AsymptoticMarcusHushChidsey(70000, 0.3)
nikms = [bv, m, amhc]

# rate constant benchmarks
# TODO: add vector of models cases
SUITE["rate constants"] = BenchmarkGroup()
SUITE["rate constants"]["non-integral models"] = BenchmarkGroup()
SUITE["rate constants"]["non-integral models"]["scalar V"] = BenchmarkGroup()
SUITE["rate constants"]["non-integral models"]["vector V"] = BenchmarkGroup()

for km in ni_kms
    modeltype = typeof(km)
    SUITE["rate constants"]["non-integral models"]["scalar V"][modeltype] = @benchmarkable compute_k(0.1, $km)
    SUITE["rate constants"]["non-integral models"]["vector V"][modeltype] = @benchmarkable compute_k(0.01:0.01:0.2, $km)
end

# TODO: this
SUITE["rate constants"]["integral models"] = BenchmarkGroup()
SUITE["rate constants"]["integral models"]["scalar V"] = BenchmarkGroup()
SUITE["rate constants"]["integral models"]["vector V"] = BenchmarkGroup()

# fit_overpotential benchmarks
# TODO: add vector cases and integral models
k = 2e3

SUITE["overpotential fitting"] = BenchmarkGroup()
SUITE["overpotential fitting"]["non-AD"] = BenchmarkGroup()
SUITE["overpotential fitting"]["AD"] = BenchmarkGroup()
for km in ni_kms
    modeltype = typeof(km)
    SUITE["overpotential fitting"]["non-AD"][modeltype] = @benchmarkable fit_overpotential($km, $k, autodiff=false)
    SUITE["overpotential fitting"]["AD"][modeltype] = @benchmarkable fit_overpotential($km, $k, autodiff=true)
end
