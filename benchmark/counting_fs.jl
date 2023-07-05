using ElectrochemicalKinetics
using TimerOutputs

bv = ButlerVolmer(300, 0.5)
mhcd_Cu = MarcusHushChidseyDOS(10000, 0.3, string("../data/DOSes/Cu_111_dos.txt"))

bv(0.01)
mhcd(0.01)
overpotential(100, bv)
overpotential(100, mhcd)

TimerOutputs.reset_timer!(ElectrochemicalKinetics.to)

function count_evals(model, arg, ad, title)
    println(title)
    overpotential(arg, model, autodiff=ad)
    show(ElectrochemicalKinetics.to)
    TimerOutputs.reset_timer!(ElectrochemicalKinetics.to);
end

# upshot of below (and a bunch of other experiments): I haven't found a case where the optimization trajectory actually diverges substantially (i.e. more than at like third/fourth sigfig) so the number of steps stays the same, and function calls is exactly half in AD case (I guess from finite differencing?)
count_evals(bv, 500, false, "BV, no autodiff")
count_evals(bv, 500, true, "BV, with autodiff")
count_evals(mhcd, 500, false, "MHCD, no autodiff")
count_evals(mhcd, 500, true, "MHCD, with autodiff")
