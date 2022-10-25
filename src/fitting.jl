using Zygote
using Optim
using NLsolve
using LinearAlgebra

# sum of squares loss in logarithmic coordinates
log_loss(y, y_pred) = (log.(y) .- log.(y_pred)).^2

# sum of squares loss in linear coordinates
linear_loss(y, y_pred) = (y .- y_pred).^2

# one that works for both signs
function janky_log_loss(y, y_pred)
    if all(y .* y_pred .> 0)
        return log_loss(abs.(y), abs.(y_pred))
    else
        # @warn "Values have different signs...defaulting to squared error"
        # this should be fine to (hopefully) push things into the right region...
        return linear_loss(y, y_pred)
    end
end

# helper fcn to make sure starting guess for `overpotential` has correct sign and to make things into nicely zippable vector inputs
function _get_guess(init_guess, k, model, a_r, a_o)
    guess = init_guess
    max_l = max(length(k), length(model), length(a_r), length(a_o))
    local k_check, a_r_check, a_o_check
    if init_guess isa Real
        guess = fill(init_guess, max_l)
    else
        if length(init_guess) != max_l
            @warn "Length mismatch between guess and required solution length"
        end
    end
    # make sure guess has correct sign
    if k isa Real
        if k == 0
            guess .= 0
        elseif k < 0
            guess[guess .> 0] .= -1 * guess[guess .> 0]
        end
        k_check = repeat([k], max_l)
    elseif k isa Vector
        guess[k .== 0] .= 0
        guess[k .< 0 .&& guess .> 0] .= -1 * guess[k .< 0 .&& guess .> 0]
        k_check = k
    end
    # NB: this assumes a_r and a_o are of the same type
    if a_r isa Real
        a_r_check = repeat([a_r], max_l)
        a_o_check = repeat([a_o], max_l)
    elseif a_r isa Vector
        a_r_check = a_r
        a_o_check = a_o
        @assert length(a_r) == length(a_o)
    end
    @assert length(a_r_check) == length(k_check) == length(guess)
    return guess, k_check, a_r_check, a_o_check
end

# TODO: support reaction direction flag here
"""
    overpotential(k, model; kwargs...)

Given values for current/rate constant and specified model parameters, find the overpotential that would have resulted in it. (This is the inverse of the `rate_constant` function.)

NOTE that this currently only solves for net reaction rates.
"""
function overpotential(k, model::KineticModel; a_r=1.0, a_o=1.0, guess = 0.1, T = 298, loss = janky_log_loss, autodiff = true, verbose=false, warn=true, lin_thresh=0.05, kwargs...)
    # wherever magnitude of k is small (less than lin_thresh * thermal voltage) we don't need to do the full solve...
    guess, k_check, a_r_check, a_o_check = _get_guess(guess, k, model, a_r, a_o)
    guess_solve = guess
    k_solve = k_check

    inds_linearize = findall(abs.(k_check) .<= lin_thresh .* rate_constant(kB*T, model; a_r=a_r_check, a_o=a_o_check, T=T, kwargs...))

    inds_solve = setdiff(1:length(k_check), inds_linearize)
    a_r_solve = a_r_check
    a_o_solve = a_o_check
    if length(inds_linearize) > 0
        k_solve = k_solve[inds_solve]
        guess_solve = guess[inds_solve]
        a_r_solve = a_r_check[inds_solve]
        a_o_solve = a_o_check[inds_solve]
    end

    sol_return = zero(guess)

    if length(k_solve) > 0
        function compare_k!(storage, V)
            storage .= loss(k, rate_constant(V, model; a_r=a_r_solve, a_o=a_o_solve, T = T, kwargs...))
        end
        Vs = if autodiff
            function myfun!(F, J, V)
                # println("V=",V)
                # println("loss=", loss(k, rate_constant(V, model; T=T, kwargs...)))
                if isnothing(J)
                    F .= loss(k_solve, rate_constant(V, model; a_r=a_r_solve, a_o=a_o_solve, T = T, kwargs...))
                elseif isnothing(F) && !isnothing(J)
                    gs = Zygote.gradient(V) do V
                        Zygote.forwarddiff(V) do V
                            loss(k_solve, rate_constant(V, model; a_r=a_r_solve, a_o=a_o_solve, T = T, kwargs...)) |> sum
                        end
                    end[1]
                    J .= diagm(gs)
                else
                    y, back = Zygote.pullback(V) do V
                        Zygote.forwarddiff(V) do V
                            loss(k_solve, rate_constant(V, model; a_r=a_r_solve, a_o=a_o_solve, T = T, kwargs...)) |> sum
                        end
                    end
                    F .= y
                    J .= diagm(back(one.(y))[1])
                end
                # println("J=", J)
            end
            Vs = nlsolve(only_fj!(myfun!), guess_solve, show_trace=verbose, extended_trace=true)
        else
            Vs = nlsolve(compare_k!, guess_solve, show_trace=verbose)
        end
        if !converged(Vs) && warn
            @warn "Overpotential fit not fully converged for $k...you may have fed in an unreachable reaction rate!"
        end
        sol = Vs.zero
        sol_return[inds_solve] = sol
    end

    # did we skip any solves? If so, we need to actually do the linearized solve and reconstruct the output array
    if length(inds_linearize) > 0
        k_lin = k_check[inds_linearize]
        slope = Zygote.jacobian(0) do V
            Zygote.forwarddiff(V) do V
                rate_constant(V, model; a_o=a_o_check[inds_linearize], a_r=a_r_check[inds_linearize], T=T, kwargs...)
            end
        end[1]
        sol_return[inds_linearize] = k_lin ./ slope
    end
    if length(sol_return) == 1
        sol_return = sol_return[1]
    end
    sol_return
end

inv!(x) = x .= inv.(x)

# TODO: test this
Zygote.@adjoint function overpotential(k, model, guess; loss = log_loss, T = 298, autodiff = true, verbose = false, kw...)
  Vs = overpotential(k, model, guess; loss=loss, verbose=verbose, autodiff=autodiff, T=T, kw...)
  function back(vs)
    gs = Zygote.jacobian(vs) do V
      Zygote.forwarddiff(V) do V
        rate_constant(V, model; T=T, kw...)
      end
    end[1]
    inv.(gs)
  end
  Vs, Δ -> (Δ .* back(Vs), nothing, nothing)
end

# basically the gradient of overpotential should just be the inverse of the gradient of the forward solve at that point, i.e. if rate_constant(V, model) ==k then overpotential(model, k)==V , so we don’t need to diff through the actual solve in overpotential because the gradient should just be the inverse of the gradient of rate_constant at V


# multiple models, one k value (used in thermo example)
# overpotential(k::Real, models::Vector{<:KineticModel}, guess=fill(0.1, length(k)); kwargs...) = overpotential.(Ref(k), models, Ref(guess); kwargs...)

fitting_params(t::Type{<:NonIntegralModel}) = fieldnames(t)
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
* `MarcusHushChidsey`: average_dos OR dos::DOSData OR dos_file::String
* `MarcusHushChidseyDOS`: dos::DOSData OR dos_file
Some are always options...
* `param_bounds::Dict{Symbol,Any}`: ranges of guesses for relevant model parameters. (must include all necessary keys, but defaults to some sensible ranges if not provided, see `default_param_bounds`...note that you should provide this for faster fitting if you know bounds)
* E_min and E_max for integral models...defaults to +/- 100kT or in case of MarcusHushChidseyDOS, to energy bounds on DOS data
"""
function fit_model(
    exp_data::Matrix,
    model_type::Type{<:KineticModel};
    param_bounds::Dict = default_param_bounds,
    T = 298,
    kwargs...
)
    V_vals = exp_data[:, 1]
    eval_model(model) = [rate_constant(V, model; T=T) for V in V_vals]
    _fit_model(exp_data, model_type, param_bounds, eval_model)
end

# TODO: check if this works with Cq
function fit_model(
    exp_data::Matrix,
    model_type::Type{<:IntegralModel};
    param_bounds::Dict = default_param_bounds,
    T = 298,
    E_min = -100 * kB*T,
    E_max = 100 * kB*T,
    kwargs...,
)
    V_vals = exp_data[:, 1]
    eval_model(model) = rate_constant(V_vals,
                                  model;
                                  T = T,
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
    sq_error(I_pred) = sum((abs.(I_pred) .- I_vals) .^ 2)

    # WTF, why do I need these here again...just when I think I understand scope
    fitting_params(t::Type{<:KineticModel}) = fieldnames(t)
    fitting_params(::Type{MarcusHushChidsey}) = (:A, :λ)
    fitting_params(::Type{MarcusHushChidseyDOS}) = (:A, :λ)

    # find best-fitting params
    opt_func = params -> sq_error(model_evaluator(model_builder(params)))
    local best_params
    # function grad!(s, x)
    #     gs = gradient(params -> opt_func(params), x)[1]
    #     for i in 1:length(x)
    #         s[i] = gs[i]
    #     end
    # end
    lower = [param_bounds[p][1] for p in fitting_params(model_type)]
    upper = [param_bounds[p][2] for p in fitting_params(model_type)]
    init_guess = 0.5 .* (lower .+ upper)

    # set optimisers based on the model type
    # IntegralModels performed better with GradientDescent
    optimizer, opts = if model_type <: IntegralModel
      opts = Optim.Options(show_trace = false,
                           iterations = 5,
                           outer_iterations = 7,
                           x_tol = 1.,
                           f_tol = 1e2)
      Fminbox(NelderMead()), opts
    else
      Fminbox(NelderMead()), Optim.Options()
    end
    opt = optimize(opt_func, lower, upper, init_guess, optimizer, opts)
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
