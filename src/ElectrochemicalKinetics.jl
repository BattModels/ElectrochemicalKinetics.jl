module ElectrochemicalKinetics

include("dos.jl")
export DOS
using .DOS: DOSData, get_dos
export DOSData, get_dos

include("kinetic_models.jl")
export ButlerVolmer, AsymptoticMarcusHushChidsey
export fermi_dirac
export integrand
export KineticModel, IntegralModel
export Marcus, MarcusHushChidsey, MarcusHushChidseyDOS

include("quantum_capacitance.jl")

include("rate_constant.jl")
export compute_k, compute_k_cq

include("fitting.jl")
include("integration_utils.jl")
export fitting_params, fit_model, fit_overpotential, is_dosmodel

include("phase_diagrams.jl")
export h, s, g_thermo, µ_thermo, µ_kinetic, g_kinetic, common_tangent, find_phase_boundaries

include("plot_fcns.jl")
export plot_models, plot_exp_and_models

end
