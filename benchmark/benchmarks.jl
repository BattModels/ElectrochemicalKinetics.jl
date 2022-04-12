using BenchmarkTools
using ElectrochemicalKinetics

const SUITE = BenchmarkGroup()

SUITE["integrals"] = BenchmarkGroup()

SUITE["integrals"]["scale"] = @benchmarkable ElectrochemicalKinetics.scale(0.1, 0.2)

# reusing the same ones from phase_diagram_tests
bv = ButlerVolmer(300)
m = Marcus(5000, 0.3)
amhc = AsymptoticMarcusHushChidsey(70000, 0.3)
kms = [bv, m, amhc]

# rate constant benchmarks
# TODO: add vector cases and integral models
SUITE["rate constants"] = BenchmarkGroup()
SUITE["rate constants"]["compute_k"] = BenchmarkGroup()
SUITE["rate cosntants"]["callable"] = BenchmarkGroup()

# fit_overpotential benchmarks
# TODO: add vector cases and integral models
k = 2e3

SUITE["overpotential fitting"] = BenchmarkGroup()
SUITE["overpotential fitting"]["non-AD"] = BenchmarkGroup()
SUITE["overpotential fitting"]["AD"] = BenchmarkGroup()
for km in kms
    modeltype = typeof(km)
    SUITE["overpotential fitting"]["non-AD"][modeltype] = @benchmarkable fit_overpotential($bv, $k, autodiff=false)
    SUITE["overpotential fitting"]["AD"][modeltype] = @benchmarkable fit_overpotential($bv, $k, autodiff=true)
end
