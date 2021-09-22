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
export fit_params, fit_model

include("plot_fcns.jl")
export plot_models, plot_exp_and_models

end
