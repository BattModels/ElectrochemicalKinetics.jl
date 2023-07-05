using Zygote
using Enzyme
using ElectrochemicalKinetics

bv = ButlerVolmer(300, 0.5)
mhcd = MarcusHushChidseyDOS(10000, 0.3, string(dirname(pathof(ElectrochemicalKinetics)), "/../data/DOSes/Cu_111_dos.txt"))

# first let's just compare time to get value and derivative for rate constants...this is all for SCALARS
function rc_zyg_f(model, V)
    y, back = Zygote.pullback(V) do V
        Zygote.forwarddiff(V) do V
            rate_constant(V, model)
        end
    end
    y, back(one.(y))[1]
end

function rc_zyg_r(model, V)
    y, back = Zygote.pullback(V) do V
        rate_constant(V, model)
    end
    y, back(one.(y))[1]
end

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
