using Zygote
using Enzyme
Enzyme.API.runtimeActivity!(true)
using BenchmarkTools
using ElectrochemicalKinetics

bv = ButlerVolmer(300, 0.5)
mhcd = MarcusHushChidseyDOS(10000, 0.3, string(dirname(pathof(ElectrochemicalKinetics)), "/../data/DOSes/Cu_111_dos.txt"))

V_test = 0.2
bv(V_test)
mhcd(V_test)

# first let's just compare time to get value and derivative for rate constants...this is all for SCALARS
function rc_zyg_f(model, V)
    y, back = Zygote.pullback(V) do V
        Zygote.forwarddiff(V) do V
            rate_constant(V, model)
        end
    end
    y, back(one.(y))[1]
end

# this one is so slow, not bothering, but keeping here for reference
# function rc_zyg_r(model, V)
#     y, back = Zygote.pullback(V) do V
#         rate_constant(V, model)
#     end
#     y, back(one.(y))[1]
# end

function rc_enz_f(model, V)
    f = V->rate_constant(V, model)
    # I think the 1.0 just multiplies the resulting gradient (for chain-rule accumulation)
    autodiff(Forward, f, Duplicated, Duplicated(V, 1.0))
end

function rc_enz_r(model, V)
    f = V->rate_constant(V, model)
    # I think the 1.0 just multiplies the resulting gradient (for chain-rule accumulation)
    deriv, primal = autodiff(ReverseWithPrimal, f, Active, Active(V))
    primal, deriv[1]
end

# spelling out a for loop here because BenchmarkTools got mad; this just has to be pasted into the REPL lol

model = bv
println(model)
println("Zygote, forward")
@benchmark rc_zyg_f($model, $V_test)
println("Enzyme, forward")
@benchmark rc_enz_f($model, $V_test)
println("Enzyme, reverse")
@benchmark rc_enz_r($model, $V_test)

model = mhcd
println(model)
println("Zygote, forward")
@benchmark rc_zyg_f($model, $V_test)
println("Enzyme, forward")
@benchmark rc_enz_f($model, $V_test)

# # this one still crashes...
# println("Enzyme, reverse")
# @benchmark rc_enz_r($model, $V_test)
