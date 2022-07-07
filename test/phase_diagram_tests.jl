using QuadGK

muoA = 0.02
muoB = 0.03
T = 298

@testset "Basic free energy functions" begin
    # enthalpy of mixing
    @test h(0.5) == 0.025 # value
    @test h(0.1) ≈ h(0.9) # symmetry
    @test h([0.1, 0.2, 0.7]) ≈ [0.009, 0.016, 0.021] # vector input
    @test h(0.5; Ω=1) == 0.25 # kwargs

    # similar set for entropy...
    @test isapprox(s(0.5), -kB * log(0.5), atol=1e-5)
    @test s(0.1) == s(0.9)
    @test s([0.1, 0.2, 0.7]) ≈ [2.80134626433e-5, 4.3121323932e-5, 5.2640192129e-5]

    # ...and thermodynamic free energy...
    @test g_thermo(0.5) == h(0.5) - T * s(0.5) + muoA * 0.5 + muoB * 0.5 ≈ 0.032200226968
    @test isapprox(g_thermo(0.0), muoA, atol=1e-8)
    @test isapprox(g_thermo(1.0), muoB, atol=1e-8)
    @test g_thermo([0.1, 0.2, 0.7]) ≈ [0.021651988, 0.02514984585595, 0.0323132232]

    # ...and thermodynamic chemical potential.
    @test all(isinf.(μ_thermo([0.0, 1.0])))
    @test μ_thermo(0.5) == μ_thermo(0.5, T=500) ≈ muoB - muoA
    @test μ_thermo(0.1) ≈  0.0335760350387
    @test μ_thermo(0.9, T=400) ≈ 0.0057368657199
    @test μ_thermo([0.1, 0.2, 0.7]) ≈ [0.0335760350387, 0.034400441691438, -0.00824168486]
    @test μ_thermo([0.1, 0.2, 0.7], T=350) ≈ [0.023730242495, 0.0281884382282, -0.004444931883]
end

# prefactors are set so that k vs. η plots roughly line up for small η
bv = ButlerVolmer(300)
m = Marcus(5000, 0.3)
amhc = AsymptoticMarcusHushChidsey(70000, 0.3)
kms = [bv, m, amhc]

xs = [0.1, 0.5, 0.95]

@testset "Kinetic μ" begin
    # all these numbers are just references evaluated as of 2/24/22
    μ_50_vals = Dict(
        :ButlerVolmer=>Dict(
            0.1 => 0.03838429151, 
            0.5 => 0.01862765755, 
            0.95 => 0.062373147), 
        :Marcus=>Dict(
            0.1 => 0.03874237118, 
            0.5 => 0.01928206901, 
            0.95 =>  0.0749319471),
        :AsymptoticMarcusHushChidsey=>Dict(
            0.1 => 0.03890969477, 
            0.5 => 0.01957132911, 
            0.95 => 0.07351263138)
        )
    μ_200_vals = Dict(
        :ButlerVolmer=>[0.05242161033, 0.04250973994, 0.13059172875], 
        :Marcus=>[0.0540075165676, 0.0458811724237, 0.212208797161],
        :AsymptoticMarcusHushChidsey=> [0.0545105175689, 0.04634317606989, 0.174525852191]
        )
    μ_100_T400_vals = Dict(
        :ButlerVolmer=>[0.0238392619734, 0.027028750035, 0.121278859754] , 
        :Marcus=>[0.0245705001602, 0.02842981631828, 0.1529969664115],
        :AsymptoticMarcusHushChidsey=>[0.02488745742498, 0.02890884116453, 0.14468534057]
        )
    for km in kms
        @testset "$(typeof(km))" begin
            μ_50 = μ_kinetic(50, km)
            μ_200 = μ_kinetic(200, km)
            μ_100_T400 = μ_kinetic(100, km, T=400)
            # test that the right (approximate) relationships hold
            μs = μ_50_vals[typeof(km).name.name]
            for x in xs
                @test μ_50(x) ≈ μs[x] ≈ μ_thermo(x) + overpotential(50, (1-x)*km)
            end

            # test vector inputs
            @test μ_200(xs) ≈ μ_200_vals[typeof(km).name.name]
            @test μ_100_T400(xs) == μ_100_T400.(xs) ≈ μ_100_T400_vals[typeof(km).name.name] 
        end
    end
end

@testset "Kinetic g" begin
    g_200_T400_vals = Dict(
        :ButlerVolmer => [0.0205857425, 0.037693617, 0.06661696], 
        :Marcus => [0.02073470111, 0.03874647481, 0.07399607032],
        :AsymptoticMarcusHushChidsey => [0.0207838185, 0.0390031065, 0.0729910686]
        )
    g_50_2_vals = Dict(:ButlerVolmer=>0.03519728018258, :Marcus=>0.03542169346909, :AsymptoticMarcusHushChidsey=>0.03552448666)
    for km in kms
        @testset "$(typeof(km))" begin
            μ_50 = μ_kinetic(50, km)
            g_50 = g_kinetic(50, km)
            g_200_T400 = g_kinetic(200, km, T=400)

            # integral of the derivative should be the original fcn (up to a constant, which we know)
            integrated_μs = [v[1] for v in quadgk.(μ_50, zero(xs), xs)]
            @test all(isapprox.(g_50(xs), integrated_μs .+ muoA, atol=1e-5))

            # check scalar input works
            @test g_50(xs[2]) == g_50(xs)[2] 
            @test isapprox(g_50(xs)[2], g_50_2_vals[typeof(km).name.name], atol=1e-5)

            # check a few other values
            @test all(isapprox.(g_200_T400(xs), g_200_T400_vals[typeof(km).name.name], atol=1e-5))
        end
    end
end

@testset "Phase Diagram" begin
    km = bv # TODO: expand this

    # simplest case, just one pair of x values (this function is still pretty slow though)
    v1 = find_phase_boundaries(100, km)

    @test all(isapprox.(common_tangent(v1, 100, km), Ref(0.0), atol=1e-6))
    v2 = find_phase_boundaries(100, km, T=350)     
    @test all(isapprox.(common_tangent(v2, 100, km, T=350), Ref(0.0), atol=1e-5))
    # they should get "narrower" with temperature
    @test v2[1] > v1[1]
    @test v2[2] < v1[2]
    v3 = find_phase_boundaries(400, km)
    # ...and also with current
    @test v3[1] > v1[1]
    @test v3[2] < v1[2]

    # test actual numerical values too
    @test isapprox(v1, [0.0458468, 0.839934], atol=1e-5)
    @test isapprox(v2, [0.08379, 0.795023], atol=1e-5)
    @test isapprox(v3, [0.1076, 0.675928], atol=1e-5)

    # TODO: next for multiple currents at once (may require some syntax tweaks)

end
