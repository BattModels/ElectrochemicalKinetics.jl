module DOS

using Interpolations
using DelimitedFiles
using Statistics

export DOSData, get_dos

"""
    DOSData(dos_file; Ef=0, cut_energy=false)

Struct for storing density of states (DOS) information. Can be constructed from a CSV file.

# Fields
* interp_func: function that returns interpolated values of DOS for any energy between E_min and E_max
* average_value: Average value of DOS over this interval
* E_min, E_max: bounds over which this DOS is defined
"""
struct DOSData
    interp_func
    average_value
    E_min
    E_max
end

DOSData(dos_file; Ef=0, cut_energy=false) = DOSData(get_dos(dos_file; Ef, cut_energy)...)

# pretty printing
Base.show(io::IO, dd::DOSData) = print(io, "DOSData: avg value $(round(dd.average_value, sigdigits=3)) from energy $(round(dd.E_min, sigdigits=2)) to $(round(dd.E_max, sigdigits=2))")

(dd::DOSData)(E::Real) = dd.interp_func(E)

import Base.length
length(::DOSData) = 1

"""
    get_dos(dos_file; Ef=0, cut_energy=false)

Retrieve density of states data from a file at path `dos_file`, assuming that the Fermi level is at energy `Ef`.

Returns a callable interpolated DOS, the average value of DOS, and the lower and upper bound energies.

# Notes
* `dos_file` is assumed to contain two colums, energy and DOS value
* If `cut_energy==true`, return data in a symmetric range of energies around `Ef` (this is sometimes useful for nice-looking plots). Note that the code presumes the cutoff for this comes from the upper bound. This would not be hard to fix, I just have not done so as yet.
"""
function get_dos(dos_file::String; kwargs...)
    # dos_data = Zygote.ignore() do
    #     x = readdlm(dos_file, Float32, skipstart=1)
    #     x
    # end
    dos_data = readdlm(dos_file, skipstart=1)
    get_dos(dos_data; kwargs...)
end
function get_dos(dd::Matrix; Ef=0, cut_energy=false)
    # dos_data = Zygote.ignore() do 
    #     x = dd
    #     # recenter so Ef=0
    #     x[:,1] = x[:,1] .- Ef
    #     x
    # end

    x = dd
    # recenter so Ef=0
    x[:,1] = x[:,1] .- Ef
    
    E_max = x[end, 1]
    if cut_energy #TODO: fix this so that it works if magnitude of upper bound is smaller
        # for now, we'll pick a symmetric energy range about 0 (i.e. the new Ef)
        len_keep = sum(x[:,1] .> -E_max)
        # cut off
        x = x[end-len_keep:end,:]
    end
    E_min = x[1,1]
    average_dos = mean(x[:,2]) # for whole structure
    # interpolate
    E_step = mean(x[2:end,1].-x[1:end-1,1])
    dos_interp = scale(interpolate(x[:,2], BSpline(Linear())), range(E_min, E_max+0.0001, step=E_step))
    return dos_interp, average_dos, E_min, E_max
end

end
