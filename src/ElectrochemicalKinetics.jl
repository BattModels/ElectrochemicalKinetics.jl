module ElectrochemicalKinetics

const kB = 8.61733326e-5

include("utils/dos.jl")
export DOS
using .DOS: DOSData, get_dos
export DOSData, get_dos

include("kinetic_models/kinetic_models.jl")
export ButlerVolmer, AsymptoticMarcusHushChidsey, LinearizedKineticModel
export fermi_dirac
export integrand
export KineticModel, NonIntegralModel, IntegralModel
export Marcus, MarcusHushChidsey, MarcusHushChidseyDOS
export rate_constant, rate_constant_cq

include("fitting.jl")
include("utils/integration_utils.jl")
export fitting_params, fit_model, overpotential, is_dosmodel

include("phase_diagrams.jl")
export kB, h, s, g_thermo, µ_thermo 
export µ_kinetic, g_kinetic
export common_tangent, find_phase_boundaries, phase_diagram

include("plot_fcns.jl")
export plot_models, plot_exp_and_models

end
