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
    @test s(0.5) ≈ -kB * 2 * 0.5 * log(0.5)
    @test s(0.1) == s(0.9)
    @test s([0.1, 0.2, 0.7]) ≈ [2.8012399817e-5, 4.3119676836e-5, 5.2638176908e-5]

    # ...and thermodynamic free energy...
    @test g_thermo(0.5) == h(0.5) - T * s(0.5) + muoA * 0.5 + muoB * 0.5 ≈ 0.032200909
    @test g_thermo(0.0) ≈ muoA
    @test g_thermo(1.0) ≈ muoB
    @test g_thermo([0.1, 0.2, 0.7]) ≈ [0.02165230485, 0.0251503363, 0.032313823]

    # ...and thermodynamic chemical potential.
    @test all(isinf.(μ_thermo([0.0, 1.0])))
    @test μ_thermo(0.5) == μ_thermo(0.5, T=500) ≈ muoB - muoA
    @test μ_thermo(0.1) ≈ 0.033578217
    @test μ_thermo(0.9, T=400) ≈ 0.0057339367
    @test μ_thermo([0.1, 0.2, 0.7]) ≈ [0.033578217, 0.034401818, -0.008242526]
    @test μ_thermo([0.1, 0.2, 0.7], T=350) ≈ [0.023732805, 0.028190055, -0.00444592]
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
        ButlerVolmer=>Dict(
            0.1 => 0.0383864736, 
            0.5 => 0.01862765755, 
            0.95 => 0.06237022288), 
        Marcus=>Dict(
            0.1 => 0.02841242588577, 
            0.5 => 0.01928206795, 
            0.95 => 0.07492901629449),
        AsymptoticMarcusHushChidsey=>Dict(
            0.1 => 0.0389118763, 
            0.5 => 0.0195713281175, 
            0.95 => 0.0735097011)
        )
    μ_200_vals = Dict(
        ButlerVolmer=>[0.0524237924, 0.04250974, 0.13058880458], 
        Marcus=>[0.054009696387679, 0.045881168639779, 0.212205862750],
        AsymptoticMarcusHushChidsey=>[0.05451269755, 0.046343172588, 0.17452291988558]
        )
    μ_100_T400_vals = Dict(
        ButlerVolmer=>[0.02384219096, 0.02702875, 0.1212749347], 
        Marcus=>[0.024573427974, 0.028429814254, 0.1529930329757],
        AsymptoticMarcusHushChidsey=>[0.024890385314, 0.028908839244, 0.14468140815]
        )
    for km in kms
        @testset "$(typeof(km))" begin
            μ_50 = μ_kinetic(50, km)
            μ_200 = μ_kinetic(200, km)
            μ_100_T400 = μ_kinetic(100, km, T=400)
            # test that the right (approximate) relationships hold
            μs = μ_50_vals[typeof(km)]
            for x in xs
                @test μ_50(x) ≈ μs[x] ≈ μ_thermo(x) + fit_overpotential(50, (1-x)*km)
            end

            # test vector inputs
            @test μ_200(xs) ≈ μ_200_vals[typeof(km)]
            @test μ_100_T400(xs) == μ_100_T400.(xs) ≈ μ_100_T400_vals[typeof(km)] 
        end
    end
end

@testset "Kinetic g" begin
    g_200_T400_vals = Dict(
        ButlerVolmer => [0.0205857425, 0.037693617, 0.06661696], 
        Marcus => [0.02073470111, 0.03874647481, 0.07399607032],
        AsymptoticMarcusHushChidsey => [0.0207838185, 0.0390031065, 0.0729910686]
        )
    g_50_2_vals = Dict(ButlerVolmer=>0.03519728018258, Marcus=>0.03395267141265, AsymptoticMarcusHushChidsey=>0.0347968044)
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
            @test isapprox(g_50(xs)[2], g_50_2_vals[typeof(km)], atol=1e-5)

            # check a few other values
            @test all(isapprox.(g_200_T400(xs), g_200_T400_vals[typeof(km)], atol=1e-5))
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
    @test isapprox(v1, [0.04584 0.839942], atol=1e-6)
    @test isapprox(v2, [0.0837787 0.795035], atol=1e-6)
    @test isapprox(v3, [0.107585 0.675948], atol=1e-6)

    # TODO: next for multiple currents at once (may require some syntax tweaks)

end
