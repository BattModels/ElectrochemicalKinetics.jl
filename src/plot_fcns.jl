using Plots
using ElectrochemicalKinetics

"""
    plot_models(models::Vector{<:KineticModel}; kwargs...)

Plot predicted rate constants for each model in the provided list.

# Keyword Arguments
* V_min and V_max: bounds of voltage, defaults to +/- 1.0
* plot_title
* kwargs for rate_constant function
"""
function plot_models(
    models::Vector{<:KineticModel};
    plot_title = "",
    V_min = -1.0,
    V_max = 1.0,
    kwargs...,
)
    V_range = collect(range(V_min, V_max, length = 200))
    xs = repeat([V_range], length(models))
    ys = map(model -> abs.(rate_constant(V_range, model; kwargs...)), models)
    plot(
        xs,
        ys,
        seriestype = [:line repeat([:line], length(models) - 1)...],
        label = [repr(models[1]) repr.(models[2:end])...],
        xlabel = "V",
        ylabel = "log(k or I)",
        yscale = :log10,
        leg = :bottomright,
        title = plot_title,
    )
end

"""
    plot_exp_and_models(exp_data::Matrix, models::Vector{<:KineticModel}; kwargs...)

Plot predicted rate constants for each model in the provided list.

# Keyword Arguments
* V_min and V_max: bounds of voltage, defaults to +/- 1.0
* plot_title
* kwargs for rate_constant function
"""
function plot_exp_and_models(
    exp_data::Matrix,
    models::Vector{<:KineticModel};
    plot_title = "",
    kwargs...,
)
    V = exp_data[:, 1]
    V_mag = 1.1 * maximum(abs.(V))
    V_range = collect(range(-V_mag, V_mag, length = 200))
    xs = Vector[V, repeat([V_range], length(models))...]
    ys = Vector[
        exp_data[:, 2],
        map(model -> abs.(rate_constant(V_range, model; kwargs...)), models)...,
    ]

    # scatter plot of experimental data, lines for fits
    plot(
        xs,
        ys,
        seriestype = [:scatter repeat([:line], length(models))...],
        label = ["experiment" repr.(models)...],
        xlabel = "V",
        ylabel = "log(k or I)",
        yscale = :log10,
        leg = :bottomright,
        title = plot_title,
    )
end

# if you want to visualize each term in MHC+DOS integrand separately
function plot_integrand_components(
    dos_file,
    E_min,
    E_max;
    eη = 0,
    kT = 0.026,
    λ = 0.26,
    length = 500,
    plot_title = "",
)
    E_range = range(E_min + 1e-9, E_max - 1e-9, length = length)
    dos = DOSData(dos_file)
    p1 = plot(
        E_range,
        dos.(E_range),
        color = :black,
        ylabel = "electrode states",
        label = "electrons/holes",
        title = plot_title,
    )

    fd_vals_e = fermi_dirac.(E_range; kT = kT)
    fd_vals_h = 1 .- fd_vals_e
    p2 = plot(
        E_range,
        hcat(fd_vals_e, fd_vals_h),
        label = ["electrons" "holes"],
        ylabel = "occupation",
        legend = :right,
    )

    marcus_model = Marcus(λ)
    ox_dos_vals = integrand(marcus_model, eη, true; kT = kT).(E_range)
    red_dos_vals = integrand(marcus_model, eη, false; kT = kT).(E_range)
    p3 = plot(
        E_range,
        hcat(red_dos_vals, ox_dos_vals),
        label = ["reduction" "oxidation"],
        ylabel = "electrolyte states",
    )

    plot(p1, p2, p3, layout = (3, 1))
end

