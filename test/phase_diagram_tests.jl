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
bv = ButlerVolmer(300, 0.5)
m = Marcus(5000, 0.3)
amhc = AsymptoticMarcusHushChidsey(70000, 0.3)
kms = [bv, m, amhc]

xs = [0.1, 0.5, 0.95]

@testset "Kinetic μ" begin
    # all these numbers are just references evaluated as of 2/24/22
    μ_50_vals = Dict(
        ButlerVolmer=>Dict(
            0.1 => 0.038327053476, 
            0.5 => 0.0185210294275, 
            0.95 => 0.061545233344), 
        Marcus=>Dict(
            0.1 => 0.03886794468, 
            0.5 => 0.019502081602, 
            0.95 => 0.0760262375),
        AsymptoticMarcusHushChidsey=>Dict(
            0.1 => 0.03908885163, 
            0.5 => 0.0198878883417, 
            0.95 => 0.0751700053714)
        )
    μ_200_vals = Dict(
        ButlerVolmer=>[0.05219089, 0.04210798, 0.1289228143], 
        Marcus=>[0.054475345584, 0.04662840026, 0.21271749102576],
        AsymptoticMarcusHushChidsey=>[0.055178277917, 0.0474096198651, 0.17595126006354]
        )
    μ_100_T400_vals = Dict(
        ButlerVolmer=>[0.023723841579, 0.02681830128, 0.120038391227], 
        Marcus=>[0.0248176854314, 0.0288548311435, 0.153980373795],
        AsymptoticMarcusHushChidsey=>[0.0252403003, 0.029515176702, 0.14635806345]
        )
    for km in kms
        @testset "$(typeof(km))" begin
            μ_50 = x -> μ_kinetic(x, 50, km)
            μ_200 = x -> μ_kinetic(x, 200, km)
            μ_100_T400 = x -> μ_kinetic(x, 100, km, T=400)
            # test that the right (approximate) relationships hold
            μs = μ_50_vals[typeof(km)]
            for x in xs
                @test μ_50(x)[1] ≈ μs[x] ≈ μ_thermo(x) + overpotential(50, km, 1-x)[1]
            end

            # test vector inputs
            @test all(isapprox.(μ_200(xs), μ_200_vals[typeof(km)], atol=1e-5))
            test100400 = µ_100_T400(xs)
            @test isapprox(test100400, [µ[1] for µ in μ_100_T400.(xs)], atol=1e-5)
            @test isapprox(test100400, μ_100_T400_vals[typeof(km)], atol=1e-5)
        end
    end
end

@testset "Kinetic g" begin
    # these cases were computed with Simpson 3/8 at a discretization of 2e-4 on 7/1/22 on the phase_diagrams_upgrade branch
    g_200_T400_vals = Dict(
        ButlerVolmer => [0.0212, 0.0415, 0.0793], 
        Marcus => [0.0201, 0.0345, 0.0619],
        AsymptoticMarcusHushChidsey => [0.01986, 0.0331, 0.05523]
        )
    g_50_2_vals = Dict(ButlerVolmer=>0.03515769, Marcus=>0.035495824, AsymptoticMarcusHushChidsey=>0.03563225)
    for km in kms
        @testset "$(typeof(km))" begin
            μ_50 = x -> μ_kinetic(x, 50, km)
            g_50 = x -> g_kinetic(x, 50, km)
            g_200_T400 = x -> g_kinetic(x, 200, km, T=400)

            # integral of the derivative should be the original fcn (up to a constant, which we know)
            integrated_μs = [v[1][1] for v in quadgk.(μ_50, zero(xs), xs)]
            @test all(isapprox.(g_50(xs), integrated_μs .+ muoA, atol=1e-4))

            # check scalar input works
            @test isapprox(g_50(xs[2]), g_50(xs)[2], atol=1e-6)
            @test isapprox(g_50(xs)[2], g_50_2_vals[typeof(km)], atol=2e-5)

            # check a few other values
            @test all(isapprox.(g_200_T400(xs), g_200_T400_vals[typeof(km)], atol=1e-3))
        end
    end
end

@testset "Phase Diagram" begin
    km = bv # TODO: expand this

    # simplest case, just one pair of x values (this function is still pretty slow though)
    v1 = find_phase_boundaries(100, km)

    @test all(isapprox.(common_tangent(v1, 100, km), Ref(0.0), atol=1e-5))
    v2 = find_phase_boundaries(100, km, T=350)
    @test all(isapprox.(common_tangent(v2, 100, km, T=350), Ref(0.0), atol=1e-3))
    # they should get "narrower" with temperature
    @test v2[1] > v1[1]
    @test v2[2] < v1[2]
    # v3 = find_phase_boundaries(400, km)
    v3 = find_phase_boundaries(120, bv)
    # ...and also with current
    @test v3[1] > v1[1]
    @test v3[2] < v1[2]

    # test actual numerical values too
    # TODO: check these better, build whole phase diagram, etc.
    # @test isapprox(v1, [0.04584 0.839942], atol=1e-6)
    # @test isapprox(v2, [0.0837787 0.795035], atol=1e-6)
    # @test isapprox(v3, [0.107585 0.675948], atol=1e-6)
    @test isapprox(v1, [0.045, 0.8405559], atol=1e-6)
    @test isapprox(v2, [0.454996, 0.790005], atol=1e-6)
    @test isapprox(v3, [0.0486987, 0.8227893], atol=1e-6)

    # TODO: next for multiple currents at once (may require some syntax tweaks)

end
