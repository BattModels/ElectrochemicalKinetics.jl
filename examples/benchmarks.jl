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
@benchmark fit_overpotential($bv, $k, autodiff=false) # ~16µs
@benchmark fit_overpotential($bv, $k, autodiff=true) # ~49µs
@benchmark fit_overpotential($m, $k) # ~30µs
@benchmark fit_overpotential($amhc, $k) #~21µs