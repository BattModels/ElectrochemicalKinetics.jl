using BlackBoxOptim

fit_params(t::Type{<:KineticModel}) = fieldnames(t)
fit_params(::Type{MarcusHushChidsey}) = (:A, :λ)
fit_params(::Type{MarcusHushChidseyDOS}) = (:A, :λ)
const default_param_bounds = Dict(:A => (0.1, 10000), :λ => (0.01, 0.5), :α => (0.01, 0.99))

"""
    fit_model(exp_data, model_type, param_bounds, kT=.026; kwargs...)
    fit_model(exp_data, model_type, param_bounds, E_min, E_max, kT=.026; kwargs...)

# Arguments
* `exp_data::Matrix`: two columns, first with voltage values, second with current
* `model_type::Type{<:KineticModel}`
* (required only if `model_type<:IntegralModel`) integral bounds `E_min` and `E_max` (note the two different method signatures for IntegralModels vs. not)

# Keyword Arguments
Requirements differ by model type...
* `ButlerVolmer`, `AsymptoticMarcusHushChidsey`: none
* `Marcus`: E_min::Float64, E_max::Float64
* `MarcusHushChidsey`: average_dos::Float64 OR dos::DOSData OR dos_file::String
* `MarcusHushChidseyDOS`: dos::DOSData OR dos_file
Some are always options...
* `param_bounds::Dict{Symbol,Any}`: ranges of guesses for relevant model parameters. (must include all necessary keys, but defaults to some sensible ranges if not provided, see `default_param_bounds`)
* some options of the `bboptimize` function from BlackBoxOptim. Default values are: MaxSteps=10000, MinDeltaFitnessTolerance=1e-9 
"""
function fit_model(
    exp_data::Matrix,
    model_type::Type{<:KineticModel};
    param_bounds::Dict = default_param_bounds,
    kT::Real = 0.026,
    kwargs...,
)
    V_vals = exp_data[:, 1]
    eval_model(model) = [compute_k(V, model; kT = kT) for V in V_vals]
    _fit_model(exp_data, model_type, param_bounds, eval_model, kwargs...)
end

function fit_model(
    exp_data::Matrix,
    model_type::Type{<:IntegralModel},
    E_min::Real,
    E_max::Real;
    param_bounds::Dict = default_param_bounds,
    kT::Real = 0.026,
    kwargs...,
)
    V_vals = exp_data[:, 1]
    eval_model(model) = [compute_k(E_min, E_max, V, model; kT = kT) for V in V_vals]
    _fit_model(exp_data, model_type, param_bounds, eval_model; kwargs...)
end

function _fit_model(
    exp_data,
    model_type::Type{<:KineticModel},
    param_bounds,
    model_evaluator;
    MaxSteps = 100000,
    MinDeltaFitnessTolerance = 1e-12,
    kwargs...,
)
    I_vals = exp_data[:, 2]

    model_builder = _get_model_builder(model_type, param_bounds; kwargs...)
    sq_error(I_pred) = sum((I_pred .- I_vals) .^ 2)

    # WTF, why do I need these here again...just when I think I understand scope
    fit_params(t::Type{<:KineticModel}) = fieldnames(t)
    fit_params(::Type{MarcusHushChidsey}) = (:A, :λ)
    fit_params(::Type{MarcusHushChidseyDOS}) = (:A, :λ)

    # find best-fitting params
    opt_func = params -> sq_error(model_evaluator(model_builder(params)))
    ss = [param_bounds[p] for p in fit_params(model_type)]
    res = bboptimize(
        opt_func;
        SearchSpace = RectSearchSpace(ss),
        NumDimensions = length(ss),
        MaxSteps = MaxSteps,
        TraceInterval = 5.0,
        MinDeltaFitnessTolerance = MinDeltaFitnessTolerance,
    )
    best_params = best_candidate(res)

    # construct and return model
    model_builder(best_params)
end

is_dosmodel(t::Type{<:KineticModel}) = false
is_dosmodel(t::Type{<:IntegralModel}) = t == MarcusHushChidsey || t == MarcusHushChidseyDOS

function _get_model_builder(model_type, param_bounds; kwargs...)
    # check that all necessary params are provided
    @assert Set(keys(param_bounds)) >= Set(fit_params(model_type))

    # check for kwargs, build DOS object if needed, etc.
    arg_names = keys(kwargs)
    local model_builder
    if is_dosmodel(model_type)
        local dos_arg
        if :dos in arg_names
            dos_arg = :dos
        elseif :dos_file in arg_names
            dos_arg = :dos_file
        elseif :average_dos in arg_names
            @assert model_type == MarcusHushChidsey "You haven't provided a DOSData object or a file from which to build it"
            dos_arg = :average_dos
        end
        model_builder = param_tuple -> model_type(param_tuple..., kwargs[dos_arg])
    else
        model_builder = param_tuple -> model_type(param_tuple...)
    end
    return model_builder
end
