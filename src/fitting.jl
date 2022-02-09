using Zygote
using Optim
using DiffImages
using NLsolve

# sum of squares loss in logarithmic coordinates
log_loss(y, y_pred) = (log.(y) .- log.(y_pred)).^2

# sum of squares loss in linear coordinates
linear_loss(y, y_pred) = (y .- y_pred).^2

"""
    fit_overpotential(model, k; kwargs...)

Given values for current/rate constant and specified model parameters, find the overpotentials that must have resulted in it.
"""
function fit_overpotential(model::KineticModel, k, forward=true; kT=.026, loss = log_loss, kwargs...)
    # start on the correct half of the Tafel plot
    if forward
        guess = 0.1
    else
        guess = -0.1
    end
    function compare_k!(storage, V)
        storage .= loss(k, compute_k(V, model; kT=kT, kwargs...))
    end

    function grad!(storage, V)
        gs = Zygote.jacobian(V -> sum(loss(k, compute_k(V, model; kT=kT, kwargs...))), V)[1]
        storage .= gs
        nothing
    end
    Vs = nlsolve(compare_k!, grad!, repeat([guess], length(k)))
    # Vs = nlsolve(compare_k!, repeat([0.1], length(k)))
    if !converged(Vs)
        @warn "Overpotential fit not fully converged...you may have fed in an unreachable reaction rate!"
    end
    Vs.zero
end

fitting_params(t::Type{<:KineticModel}) = fieldnames(t)
fitting_params(::Type{MarcusHushChidsey}) = (:A, :λ)
fitting_params(::Type{MarcusHushChidseyDOS}) = (:A, :λ)
const default_param_bounds = Dict(:A => (0.1, 50000), :λ => (0.01, 0.5), :α => (0.01, 0.99))

"""
    fit_model(exp_data, model_type; kwargs...)

# Arguments
* `exp_data::Matrix`: two columns, first with voltage values, second with current
* `model_type::Type{<:KineticModel}`

# Keyword Arguments
Requirements differ by model type...
* `ButlerVolmer`, `AsymptoticMarcusHushChidsey`, `Marcus`: none
* `MarcusHushChidsey`: average_dos::Float64 OR dos::DOSData OR dos_file::String
* `MarcusHushChidseyDOS`: dos::DOSData OR dos_file
Some are always options...
* `param_bounds::Dict{Symbol,Any}`: ranges of guesses for relevant model parameters. (must include all necessary keys, but defaults to some sensible ranges if not provided, see `default_param_bounds`...note that you should provide this for faster fitting if you know bounds)
* E_min and E_max for integral models...defaults to +/- 100kT or in case of MarcusHushChidseyDOS, to energy bounds on DOS data
"""
function fit_model(
    exp_data::Matrix,
    model_type::Type{<:KineticModel};
    param_bounds::Dict = default_param_bounds,
    kT::Real = 0.026,
    kwargs...
)
    V_vals = exp_data[:, 1]
    eval_model(model) = [compute_k(V, model; kT = kT) for V in V_vals]
    _fit_model(exp_data, model_type, param_bounds, eval_model)
end

# TODO: check if this works with Cq
function fit_model(
    exp_data::Matrix,
    model_type::Type{<:IntegralModel};
    param_bounds::Dict = default_param_bounds,
    kT::Real = 0.026,
    E_min = -100 * kT,
    E_max = 100 * kT,
    kwargs...,
)
    V_vals = exp_data[:, 1]
    eval_model(model) = compute_k(V_vals,
                                  model;
                                  kT = kT,
                                  E_min = E_min,
                                  E_max = E_max,
                                  kwargs...)
    _fit_model(exp_data, model_type, param_bounds, eval_model; kwargs...)
end

function _fit_model(
    exp_data,
    model_type::Type{<:KineticModel},
    param_bounds,
    model_evaluator;
    kwargs...,
)
    I_vals = exp_data[:, 2]

    model_builder = _get_model_builder(model_type, param_bounds; kwargs...)
    sq_error(I_pred) = sum((I_pred .- I_vals) .^ 2)

    # WTF, why do I need these here again...just when I think I understand scope
    fitting_params(t::Type{<:KineticModel}) = fieldnames(t)
    fitting_params(::Type{MarcusHushChidsey}) = (:A, :λ)
    fitting_params(::Type{MarcusHushChidseyDOS}) = (:A, :λ)

    # find best-fitting params
    # Zygote is tripped up by QuadGK, so that one has to be done with black-box optimization, but the non-integral models work with autodiff
    opt_func = params -> sq_error(model_evaluator(model_builder(params)))
    local best_params
    function grad!(s, x)
        gs = gradient(params -> opt_func(params), x)[1]
        for i in 1:length(x)
            s[i] = gs[i]
        end
    end
    lower = Float64.([param_bounds[p][1] for p in fitting_params(model_type)])
    upper = Float64.([param_bounds[p][2] for p in fitting_params(model_type)])
    init_guess = 0.5 .* (lower .+ upper)

    # set optimisers based on the model type
    # IntegralModels performed better with GradientDescent
    optimizer, opts = if model_type <: IntegralModel
      opts = Optim.Options(show_trace = false,
                           iterations = 3,
                           outer_iterations = 4,
                           x_tol = 1.,
                           f_tol = 1e3)
      Fminbox(GradientDescent()), opts
    else
      Fminbox(LBFGS()), Optim.Options()
    end
    opt = optimize(opt_func, grad!, lower, upper, init_guess, optimizer, opts)
    best_params = Optim.minimizer(opt)

    # construct and return model
    model_builder(best_params)
end

"""
    is_dosmodel(type{<:KineticModel})

Returns true if the model needs DOS information.
"""
is_dosmodel(t::Type{<:KineticModel}) = false
is_dosmodel(t::Type{<:IntegralModel}) = t == MarcusHushChidsey || t == MarcusHushChidseyDOS

function _get_model_builder(model_type, param_bounds; kwargs...)
    # check that all necessary params are provided
    @assert Set(keys(param_bounds)) >= Set(fitting_params(model_type))

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
        else
            throw(ArgumentError("I need DOS information for this model type!"))
        end
        model_builder = param_tuple -> model_type(param_tuple..., kwargs[dos_arg])
    else
        model_builder = param_tuple -> model_type(param_tuple...)
    end
    return model_builder
end
