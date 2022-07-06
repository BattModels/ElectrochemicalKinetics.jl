using ElectrochemicalKinetics
using TimerOutputs

bv = ButlerVolmer(300, 0.5)

# run once to force compilation
find_phase_boundaries(10, bv)

reset_timer!(ElectrochemicalKinetics.to)

find_phase_boundaries(12, bv)
show(ElectrochemicalKinetics.to)