include("analysis_fcns.jl")

stepsize = 0.001
lambda=0.82
r = 0.45
lambda = 0.82
V_app = range(-r, r, step=0.001)

## bernal stacked graphene
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/bernal_graphene.txt"; Ef=0)
v, v_min, v_max = get_vh("figure_data/bernal_graphene_vh_vs_vappl.txt")
V_h = [v(V) for V in V_app]
V_q = V_app - V_h
red = []
ox = []
mhc = []
for i=1:length(V_app)
    push!(mhc, compute_k_MHC_DOS(-r, r, lambda, V_app[i], dos_f; kT=.026, vq=V_q[i]))
end
open("figure_data/extended_range_DOS/k_vs_V_0_red.txt", "w") do io
    write(io, "V, MHC_DOS_k_red\n")
    writedlm(io, hcat(V_app, red))
end;
open("figure_data/extended_range_DOS/k_vs_V_0_ox.txt", "w") do io
    write(io, "V, MHC_DOS_k_ox\n")
    writedlm(io, hcat(V_app, ox))
end
open("figure_data/extended_range_DOS/k_vs_V_0_vq.txt", "w") do io
    write(io, "V, MHC_DOS_k\n")
    writedlm(io, hcat(V_app, mhc))
end

## TBLG DOS 0.42
EF=-0.0016
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/extended_range_DOS/dos_0.42.txt"; Ef=EF)
v, v_min, v_max = get_vh("figure_data/cq_calc/vh_vs_vappl_0.42.txt")
V_h = [v(V) for V in V_app]
V_q = V_app - V_h
red = []
ox = []
mhc = []
for i=1:length(V_app)
    push!(red, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=false, kT=.026, vq=V_q[i], E0=E_0))
    push!(ox, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=true, kT=.026, vq=V_q[i], E0=E_0))
    push!(mhc, compute_k_MHC_DOS(-r, r, lambda, V_app[i], dos_f; kT=.026, vq=V_q[i], E0=E_0))
end
open("figure_data/extended_range_DOS/k_vs_V_0.42_red.txt", "w") do io
    write(io, "V, MHC_DOS_k_red\n")
    writedlm(io, hcat(V_app, red))
end;
open("figure_data/extended_range_DOS/k_vs_V_0.42_ox.txt", "w") do io
    write(io, "V, MHC_DOS_k_ox\n")
    writedlm(io, hcat(V_app, ox))
end
open("figure_data/extended_range_DOS/k_vs_V_0.42_vq.txt", "w") do io
    write(io, "V, MHC_DOS_k\n")
    writedlm(io, hcat(V_app, mhc))
end

## TBLG DOS 0.77
EF=-0.0054
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/extended_range_DOS/dos_0.77.txt"; Ef=EF)
v, v_min, v_max = get_vh("figure_data/cq_calc/vh_vs_vappl_0.77.txt")
V_h = [v(V) for V in V_app]
V_q = V_app - V_h
red = []
ox = []
mhc = []
for i=1:length(V_app)
    push!(red, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=false, kT=.026, vq=V_q[i], E0=E_0))
    push!(ox, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=true, kT=.026, vq=V_q[i],E0=E_0))
    push!(mhc, compute_k_MHC_DOS(-r, r, lambda, V_app[i], dos_f; kT=.026, vq=V_q[i], E0=E_0))
end
open("figure_data/extended_range_DOS/k_vs_V_0.77_red.txt", "w") do io
    write(io, "V, MHC_DOS_k_red\n")
    writedlm(io, hcat(V_app, red))
end;
open("figure_data/extended_range_DOS/k_vs_V_0.77_ox.txt", "w") do io
    write(io, "V, MHC_DOS_k_ox\n")
    writedlm(io, hcat(V_app, ox))
end
open("figure_data/extended_range_DOS/k_vs_V_0.77_vq.txt", "w") do io
    write(io, "V, MHC_DOS_k\n")
    writedlm(io, hcat(V_app, mhc))
end

## TBLG DOS 1.15
EF=-0.0081
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/extended_range_DOS/dos_1.15.txt"; Ef=EF)
v, v_min, v_max = get_vh("figure_data/cq_calc/vh_vs_vappl_1.15.txt")
V_h = [v(V) for V in V_app]
V_q = V_app - V_h
red = []
ox = []
mhc = []
for i=1:length(V_app)
    push!(red, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=false, kT=.026, vq=V_q[i], E0=E_0))
    push!(ox, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=true, kT=.026, vq=V_q[i], E0=E_0))
    push!(mhc, compute_k_MHC_DOS(-r, r, lambda, V_app[i], dos_f; kT=.026, vq=V_q[i], E0=E_0))
end
open("figure_data/extended_range_DOS/k_vs_V_1.15_red.txt", "w") do io
    write(io, "V, MHC_DOS_k_red\n")
    writedlm(io, hcat(V_app, red))
end;
open("figure_data/extended_range_DOS/k_vs_V_1.15_ox.txt", "w") do io
    write(io, "V, MHC_DOS_k_ox\n")
    writedlm(io, hcat(V_app, ox))
end
open("figure_data/extended_range_DOS/k_vs_V_1.15_vq.txt", "w") do io
    write(io, "V, MHC_DOS_k\n")
    writedlm(io, hcat(V_app, mhc))
end

## TBLG DOS 1.34
EF=-0.009
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/extended_range_DOS/dos_1.34.txt"; Ef=EF)
v, v_min, v_max = get_vh("figure_data/cq_calc/vh_vs_vappl_1.34.txt")
V_h = [v(V) for V in V_app]
V_q = V_app - V_h
red = []
ox = []
mhc = []
for i=1:length(V_app)
    push!(red, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=false, kT=.026, vq=V_q[i], E0=E_0))
    push!(ox, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=true, kT=.026, vq=V_q[i], E0=E_0))
    push!(mhc, compute_k_MHC_DOS(-r, r, lambda, V_app[i], dos_f; kT=.026, vq=V_q[i], E0=E_0))
end
open("figure_data/extended_range_DOS/k_vs_V_1.34_red.txt", "w") do io
    write(io, "V, MHC_DOS_k_red\n")
    writedlm(io, hcat(V_app, red))
end;
open("figure_data/extended_range_DOS/k_vs_V_1.34_ox.txt", "w") do io
    write(io, "V, MHC_DOS_k_ox\n")
    writedlm(io, hcat(V_app, ox))
end
open("figure_data/extended_range_DOS/k_vs_V_1.34_vq.txt", "w") do io
    write(io, "V, MHC_DOS_k\n")
    writedlm(io, hcat(V_app, mhc))
end

## TBLG DOS 1.5
EF=-0.0097
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/extended_range_DOS/dos_1.5.txt"; Ef=EF)
v, v_min, v_max = get_vh("figure_data/cq_calc/vh_vs_vappl_1.5.txt")
V_h = [v(V) for V in V_app]
V_q = V_app - V_h
red = []
ox = []
mhc = []
for i=1:length(V_app)
    push!(red, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=false, kT=.026, vq=V_q[i], E0=E_0))
    push!(ox, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=true, kT=.026, vq=V_q[i], E0=E_0))
    push!(mhc, compute_k_MHC_DOS(-r, r, lambda, V_app[i], dos_f; kT=.026, vq=V_q[i], E0=E_0))
end
open("figure_data/extended_range_DOS/k_vs_V_1.5_red.txt", "w") do io
    write(io, "V, MHC_DOS_k_red\n")
    writedlm(io, hcat(V_app, red))
end;
open("figure_data/extended_range_DOS/k_vs_V_1.5_ox.txt", "w") do io
    write(io, "V, MHC_DOS_k_ox\n")
    writedlm(io, hcat(V_app, ox))
end
open("figure_data/extended_range_DOS/k_vs_V_1.5_vq.txt", "w") do io
    write(io, "V, MHC_DOS_k\n")
    writedlm(io, hcat(V_app, mhc))
end

## TBLG DOS 2.39
EF=-0.0113
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/extended_range_DOS/dos_2.39.txt"; Ef=EF)
v, v_min, v_max = get_vh("figure_data/cq_calc/vh_vs_vappl_2.39.txt")
V_h = [v(V) for V in V_app]
V_q = V_app - V_h
red = []
ox = []
mhc = []
for i=1:length(V_app)
    push!(red, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=false, kT=.026, vq=V_q[i], E0=E_0))
    push!(ox, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=true, kT=.026, vq=V_q[i], E0=E_0))
    push!(mhc, compute_k_MHC_DOS(-r, r, lambda, V_app[i], dos_f; kT=.026, vq=V_q[i], E0=E_0))
end
open("figure_data/extended_range_DOS/k_vs_V_2.39_red.txt", "w") do io
    write(io, "V, MHC_DOS_k_red\n")
    writedlm(io, hcat(V_app, red))
end;
open("figure_data/extended_range_DOS/k_vs_V_2.39_ox.txt", "w") do io
    write(io, "V, MHC_DOS_k_ox\n")
    writedlm(io, hcat(V_app, ox))
end
open("figure_data/extended_range_DOS/k_vs_V_2.39_vq.txt", "w") do io
    write(io, "V, MHC_DOS_k\n")
    writedlm(io, hcat(V_app, mhc))
end

## TBLG DOS 2.6
EF=-0.0115
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/extended_range_DOS/dos_2.6.txt"; Ef=EF)
v, v_min, v_max = get_vh("figure_data/cq_calc/vh_vs_vappl_2.6.txt")
V_h = [v(V) for V in V_app]
V_q = V_app - V_h
red = []
ox = []
mhc = []
for i=1:length(V_app)
    push!(red, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=false, kT=.026, vq=V_q[i], E0=E_0))
    push!(ox, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=true, kT=.026, vq=V_q[i], E0=E_0))
    push!(mhc, compute_k_MHC_DOS(-r, r, lambda, V_app[i], dos_f; kT=.026, vq=V_q[i], E0=E_0))
end
open("figure_data/extended_range_DOS/k_vs_V_2.6_red.txt", "w") do io
    write(io, "V, MHC_DOS_k_red\n")
    writedlm(io, hcat(V_app, red))
end;
open("figure_data/extended_range_DOS/k_vs_V_2.6_ox.txt", "w") do io
    write(io, "V, MHC_DOS_k_ox\n")
    writedlm(io, hcat(V_app, ox))
end
open("figure_data/extended_range_DOS/k_vs_V_2.6_vq.txt", "w") do io
    write(io, "V, MHC_DOS_k\n")
    writedlm(io, hcat(V_app, mhc))
end

## TBLG DOS 3
EF=-0.0118
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/extended_range_DOS/dos_3.txt"; Ef=EF)
v, v_min, v_max = get_vh("figure_data/cq_calc/vh_vs_vappl_3.txt")
V_h = [v(V) for V in V_app]
V_q = V_app - V_h
red = []
ox = []
mhc = []
for i=1:length(V_app)
    push!(red, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=false, kT=.026, vq=V_q[i], E0=E_0))
    push!(ox, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=true, kT=.026, vq=V_q[i], E0=E_0))
    push!(mhc, compute_k_MHC_DOS(-r, r, lambda, V_app[i], dos_f; kT=.026, vq=V_q[i], E0=E_0))
end
open("figure_data/extended_range_DOS/k_vs_V_3_red.txt", "w") do io
    write(io, "V, MHC_DOS_k_red\n")
    writedlm(io, hcat(V_app, red))
end;
open("figure_data/extended_range_DOS/k_vs_V_3_ox.txt", "w") do io
    write(io, "V, MHC_DOS_k_ox\n")
    writedlm(io, hcat(V_app, ox))
end
open("figure_data/extended_range_DOS/k_vs_V_3_vq.txt", "w") do io
    write(io, "V, MHC_DOS_k\n")
    writedlm(io, hcat(V_app, mhc))
end

## TBLG DOS 5
EF=-0.0123
dos_f, avg_dos, E_min, E_max = get_dos("DOSes/extended_range_DOS/dos_5.txt"; Ef=EF)
v, v_min, v_max = get_vh("figure_data/cq_calc/vh_vs_vappl_5.txt")
V_h = [v(V) for V in V_app]
V_q = V_app - V_h
red = []
ox = []
mhc = []
for i=1:length(V_app)
    push!(red, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=false, kT=.026, vq=V_q[i], E0=E_0))
    push!(ox, compute_k_MHC_DOS_redox(-r, r, lambda, V_app[i], dos_f; ox=true, kT=.026, vq=V_q[i], E0=E_0))
    push!(mhc, compute_k_MHC_DOS(-r, r, lambda, V_app[i], dos_f; kT=.026, vq=V_q[i], E0=E_0))
end
open("figure_data/extended_range_DOS/k_vs_V_5_red.txt", "w") do io
    write(io, "V, MHC_DOS_k_red\n")
    writedlm(io, hcat(V_app, red))
end;
open("figure_data/extended_range_DOS/k_vs_V_5_ox.txt", "w") do io
    write(io, "V, MHC_DOS_k_ox\n")
    writedlm(io, hcat(V_app, ox))
end
open("figure_data/extended_range_DOS/k_vs_V_5_vq.txt", "w") do io
    write(io, "V, MHC_DOS_k\n")
    writedlm(io, hcat(V_app, mhc))
end

