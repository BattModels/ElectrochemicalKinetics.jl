module DOS

using Interpolations
export DOSData, get_dos

struct DOSData
    interp_func
    average_value::Float64
    E_min::Float64
    E_max::Float64
end

DOSData(dos_file; Ef=0, cut_energy=false) = DOSData(get_dos(dos_file; Ef, cut_energy)...)

"""
    get_dos(dos_file; Ef=0, cut_energy=false)

Retrieve density of states data from a file at path `dos_file`, assuming that the Fermi level is at energy `Ef`.

Returns a callable interpolated DOS, the average value of DOS, and the lower and upper bound energies.

# Notes
* `dos_file` is assumed to contain two colums, energy and DOS value
* If `cut_energy==true`, return data in a symmetric range of energies around `Ef` (this is sometimes useful for nice-looking plots). Note that the code presumes the cutoff for this comes from the upper bound. This would not be hard to fix, I just have not done so as yet.
"""
function get_dos(dos_file; Ef=0, cut_energy=false)
    dos_data = readdlm(dos_file, Float64, skipstart=1)
    # recenter so Ef=0
    dos_data[:,1] = dos_data[:,1] .- Ef
    E_max = dos_data[end,1]
    if cut_energy #TODO: fix this so that it works if magnitude of upper bound is smaller
        # for now, we'll pick a symmetric energy range about 0 (i.e. the new Ef)
        len_keep = sum(dos_data[:,1] .> -E_max)
        # cut off
        dos_data = dos_data[end-len_keep:end,:]
    end
    E_min = dos_data[1,1]
    average_dos = mean(dos_data[:,2]) # for whole structure
    # interpolate
    E_step = mean(dos_data[2:end,1].-dos_data[1:end-1,1])
    dos_interp = scale(interpolate(dos_data[:,2], BSpline(Linear())), range(E_min, E_max+0.0001, step=E_step))
    return dos_interp, average_dos, E_min, E_max
end

end