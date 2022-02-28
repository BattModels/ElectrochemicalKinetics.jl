using ElectrochemicalKinetics
using BenchmarkTools

# reusing the ones from thermo_tests
bv = ButlerVolmer(300)
m = Marcus(5000, 0.3)
amhc = AsymptoticMarcusHushChidsey(70000, 0.3)
kms = [bv, m, amhc]

### fit_overpotential
k = 2e3

# benchmarks, mean times
@benchmark fit_overpotential($bv, $k) # ~16µs without AD, 500 with
@benchmark fit_overpotential($m, $k) # ~30µs, 1.25 ms
@benchmark fit_overpotential($amhc, $k) #~21µs, 690µs

# vector of models
bvs = [bv, 0.8*bv, 1.2*bv]
@benchmark fit_overpotential($bvs, $k)