using Interpolations
using MAT
using ElectrochemicalKinetics

theta_list = [0.42, 0.77, 1.15, 1.34, 1.5, 2.39, 2.6, 3, 5];

# Interpolate Fermi level
ef_data=[[0.22 -2.7201e-18];[0.43 -0.0017];[0.57 -0.0042];[0.8 -0.0057];[1.1 -0.0078];[1.2 -0.0084];[1.36 -0.0091];[1.55 -0.0098];[1.96 -0.0108];[2.36 -0.0113];[2.67 -0.0116];[2.8 -0.0117];[3.2 -0.0119];[3.73 -0.0121];[4.4 -0.0122];[5 -0.0123]];
ef_func = LinearInterpolation(ef_data[:,1], ef_data[:,2])

# read in experimental data and drop rows for which we don't have a DOS
exp_data = read(matopen("../data/babar_data/ko_exp_paper.mat"))["ko_exp_paper"]

# "fit" prefactor to 5.0° twist data, assuming an overpotential of 0.07 V and reorganization energy of 0.82 eV
A = compute_k(0.07, MarcusHushChidseyDOS(1.0, 0.82, "../data/babar_data/dos_AB_5.0.txt")) / exp_data[12,2]

# now fit overpotential at each twist angle...with the A determined as above, this currently fails anywhere remotely near the magic angle :(
η_vals = []
for theta in theta_list
    ef = ef_func(theta)
    ldos = readdlm("../data/babar_data/dos_AB_" * string(theta) * ".txt")
    mhcd = MarcusHushChidseyDOS(A, 0.82, ldos, Ef=ef)
    append!(η_vals, fit_overpotential(mhcd, exp_data[4,2]))
end