include("functions.jl")

# semiconductor case, equilibrium
dos_file = "zeeshan_files/relax_pris.dos"
dos_f, average_dos, min_E, max_E = get_dos(dos_file; Ef=0.377)

for λ in [.26, .52, .78]
    for eη in [0.0, 1.0, -1.0, -2.0]
        plot_integrand_components(dos_f, -5, 5, average_dos; eη=eη, λ=λ, plot_title=string("eη=", eη, " λ=", λ))
        savefig(string("figs/integrand_components/lambda",λ,"_pot",eη,".png"))
    end
end
