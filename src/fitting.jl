# NOW FOR THE FITTING STUFF

# goal: minimize MSE between experimental and modeled data by choosing best-fit λ
# assumes `exp_data` is an n x 2 array of voltage, current values
# if `dos_func` flag is true, `dos` param will be assumed to be a callable
# otherwise, assumed a float (average dos)
# `p` is a two-element vector with λ plus overall scale
function sq_error(p, E_min, E_max, dos, exp_data, dos_is_func=false)
    # idiot check
    if dos_is_func
        if isempty(methods(dos)) # it's not callable...
            error("It seems you haven't passed a callable DOS. Maybe you meant to set `dos_func`=false?")
        else
            pred_func = V -> p[2]*compute_k_MHC_DOS(E_min, E_max, p[1], V, dos)
        end
    else
        if !isempty(methods(dos[1])) # it's not callable...
            error("It seems you passed a callable DOS when the function expected an average. Maybe you meant to set `dos_func`=true?")
        else
            pred_func = V -> p[2]*compute_k_MHC(E_min, E_max, p[1], V, dos)
        end
    end
    # actually calculate
    V = exp_data[:,1]
    I = exp_data[:,2]
    #sum(log.((pred_func.(V) .- I).^2))
    sum((pred_func.(V) .- I).^2)
end

function best_fit_params(exp_data, dos, E_min, E_max, dos_is_func=false; λ_bounds=(0.01, 0.5), A_bounds=(0.1, 10000))
    V_vals = exp_data[:,1]
    opt_func = p -> sq_error(p, E_min-V_vals[1], E_max-V_vals[end], dos, exp_data, dos_is_func)
    res = bboptimize(opt_func; SearchSpace=RectSearchSpace([λ_bounds, A_bounds]), NumDimensions=2, MaxSteps=100000, TraceInterval=5.0, MinDeltaFitnessTolerance=1e-12)
    best_candidate(res)
end
