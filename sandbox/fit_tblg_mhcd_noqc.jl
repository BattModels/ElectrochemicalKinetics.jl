using Interpolations
using MAT
using ElectrochemicalKinetics
using Plots
using DelimitedFiles

data_path = dirname(pathof(ElectrochemicalKinetics))*"/../data/babar_data/"
dos_prefix = "dos_" # change to "dos_AA_" or "dos_AB_" as desired

theta_list = [0.42, 0.77, 1.15, 1.34, 1.5, 2.39, 2.6, 3, 5];

# Interpolate Fermi level
ef_data=[[0.22 -2.7201e-18];[0.43 -0.0017];[0.57 -0.0042];[0.8 -0.0057];[1.1 -0.0078];[1.2 -0.0084];[1.36 -0.0091];[1.55 -0.0098];[1.96 -0.0108];[2.36 -0.0113];[2.67 -0.0116];[2.8 -0.0117];[3.2 -0.0119];[3.73 -0.0121];[4.4 -0.0122];[5 -0.0123]];
ef_func = LinearInterpolation(ef_data[:,1], ef_data[:,2])

# interpolate experimental data because 0.42° is missing
exp_data = read(matopen(data_path*"ko_exp_paper.mat"))["ko_exp_paper"]
k_func = LinearInterpolation(exp_data[:,1], exp_data[:,2])


# "fit" prefactor to 5.0° twist data, assuming an overpotential of 0.07 V and reorganization energy of 0.82 eV
A = exp_data[12,2] / compute_k(0.07, MarcusHushChidseyDOS(1.0, 0.82, data_path*dos_prefix*"5.0.txt"))

# now fit overpotential at each twist angle
η_vals = []
for theta in theta_list
    ldos = readdlm(data_path * dos_prefix * string(theta) * ".txt");
    mhcd = MarcusHushChidseyDOS(A, 0.82, ldos, Ef=ef_func(theta))
    append!(η_vals, fit_overpotential(mhcd, k_func(theta)))
end

plot(theta_list, η_vals)
xlabel!("Twist Angle [°]")
ylabel!("Overpotential [V]")
