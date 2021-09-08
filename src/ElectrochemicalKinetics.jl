module ElectrochemicalKinetics

include("dos.jl")
export DOS
using .DOS: DOSData, get_dos
export DOSData, get_dos

include("quantum_capacitance.jl")

include("rate_constant.jl")
export compute_k

include("fitting.jl")

end
