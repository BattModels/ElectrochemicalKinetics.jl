using MAT
using CSV
using DelimitedFiles

theta_list = [0.42, 0.77, 1.15, 1.34, 1.5, 2.39, 2.6, 3, 5];

dos_loc = "/Volumes/Data/GDrive/CMU/research/MHC_DOS/Data_24kpts/"

for theta in theta_list
    local theta_str
    if round(theta) == theta
        theta_str = string(Int(theta))
    else
        theta_str = string(theta)
    end

    file = matopen(dos_loc*"ldos-"*theta_str*"_90x90.mat")
    data = read(file, "data");

    AA_data = data[:,1,1]
    AB_data = data[:,12,12]

    d = readdlm(dos_loc*"dos-"*theta_str*".csv", ',', Float64)
    Elist = d[:,1]

    open("../data/babar_data/dos_AA_"*string(theta)*".txt", "w") do io
        writedlm(io, [Elist AA_data])
    end
    open("../data/babar_data/dos_AB_"*string(theta)*".txt", "w") do io
        writedlm(io, [Elist AB_data])
    end
end