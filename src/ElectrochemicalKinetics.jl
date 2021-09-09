module ElectrochemicalKinetics

include("dos.jl")
export DOS
using .DOS: DOSData, get_dos
export DOSData, get_dos

include("kinetic_models.jl")
export integrand
export Marcus, MarcusHushChidsey, MarcusHushChidseyDOS

include("quantum_capacitance.jl")

include("rate_constant.jl")
export compute_k

include("fitting.jl")

end
