using Plots
using ElectrochemicalKinetics

# plotting shortcut to make MHC and MHC+DOS curves with some fixed parameters
function plot_comparison(
    dos_file,
    eη_max,
    E_min,
    E_max;
    kT = 0.026,
    λ = 0.26,
    length = 500,
    plot_title = "",
)
    eη_range = range(-eη_max, eη_max, length = length)
    MHC_model = MarcusHushChidsey(λ, dos_file)
    MHC_DOS_model = MarcusHushChidseyDOS(λ, dos_file)
    MHC_k = [compute_k(E_min, E_max, eη, MHC_model; kT = kT) for eη in eη_range]
    MHC_DOS_k = [compute_k(E_min, E_max, eη, MHC_DOS_model; kT = kT) for eη in eη_range]
    plot(
        eη_range,
        hcat(log10.(abs.(MHC_k)), log10.(abs.(MHC_DOS_k))),
        label = ["MHC" "MHC+DOS"],
        legend = :bottomright,
        xlabel = "eη",
        ylabel = "log(k)",
        title = plot_title,
        color = [:purple :green],
    )
end

# if you want to visualize each term in integrand separately
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

function plot_fits(
    exp_data,
    dos_file,
    E_min,
    E_max,
    MHC_λ,
    MHC_A,
    MHC_DOS_λ,
    MHC_DOS_A;
    plot_title = "",
)
    # compute k's
    V = exp_data[:, 1]
    V_mag = 1.1 * maximum(abs.(V))
    V_range = range(-V_mag, V_mag, length = 200)
    MHC_model = MarcusHushChidsey(MHC_λ, dos_file)
    MHC_DOS_model = MarcusHushChidseyDOS(MHC_DOS_λ, dos_file)
    MHC_k = [MHC_A * compute_k(E_min, E_max, V, MHC_model) for V in V_range]
    MHC_DOS_k = [MHC_DOS_A * compute_k(E_min, E_max, V, MHC_DOS_model) for V in V_range]

    # scatter plot of experimental data, lines for fits
    xs = Vector[V, V_range, V_range]
    ys = Vector[exp_data[:, 2], MHC_k, MHC_DOS_k]
    plot(
        xs,
        ys,
        seriestype = [:scatter :line :line],
        label = ["experiment" string(
            "MHC: λ=",
            round(MHC_λ, digits = 3),
            "; A=",
            round(MHC_A, digits = 2),
        ) string(
            "MHC+DOS: λ=",
            round(MHC_DOS_λ, digits = 3),
            "; MHC_DOS_A=",
            round(MHC_DOS_A, digits = 2),
        )],
        xlabel = "V",
        ylabel = "log(k or I)",
        yscale = :log10,
        leg = :bottomright,
        title = plot_title,
    )
end
