using QuadGK

muoA = 0.02
muoB = 0.03
T = 298

# TODO: add tests where we change Ω, etc.
# TODO: add tests for deintercalation case

@testset "Basic free energy functions" begin
    # enthalpy of mixing
    @test h(0.5) == 0.025 # value
    @test h(0.1) ≈ h(0.9) # symmetry
    @test all(h([0.1, 0.2, 0.7]) .≈ [0.009, 0.016, 0.021]) # vector input
    @test h(0.5; Ω=1) == 0.25 # kwargs

    # similar set for entropy...
    @test isapprox(s(0.5f0), -kB * log(0.5f0), atol=1f-5)
    @test s(0.1) == s(0.9)
    @test all(s([0.1, 0.2, 0.7]) .≈ [2.8013f-5, 4.3121f-5, 5.264f-5])

    # ...and thermodynamic free energy...
    @test g_thermo(0.5) == h(0.5) - T * s(0.5) + muoA * 0.5 + muoB * 0.5 ≈ 0.0322002f0
    @test isapprox(g_thermo(0.0), muoA, atol=1e-8)
    @test isapprox(g_thermo(1.0), muoB, atol=1e-8)
    @test all(g_thermo([0.1, 0.2, 0.7]) .≈ [0.021651988f0, 0.02514984585595f0, 0.0323132232f0])

    # ...and thermodynamic chemical potential.
    @test all(isinf.(μ_thermo([0.0, 1.0])))
    @test μ_thermo(0.5) == μ_thermo(0.5, T=500) ≈ muoB - muoA
    @test μ_thermo(0.1) ≈ 0.033576035f0
    @test isapprox(μ_thermo(0.9, T=400), 0.00573686571, atol=5e-6)
    @test all(μ_thermo([0.1, 0.2, 0.7]) .≈ [0.033576f0, 0.03440044f0, -0.0082417f0])
    @test all(μ_thermo([0.1, 0.2, 0.7], T=350) .≈ [0.02373024f0, 0.0281884f0, -0.00444493f0])
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
        :ButlerVolmer=>Dict(
            0.1 => 0.038325055f0, 
            0.5 => 0.018521359f0, 
            0.95 =>  0.0615507f0), 
        :Marcus=>Dict(
            0.1 => 0.03886537f0, 
            0.5 =>  0.01950138f0, 
            0.95 =>  0.0760257f0),
        :AsymptoticMarcusHushChidsey=>Dict(
            0.1 => 0.0390861f0, 
            0.5 => 0.01988688f0, 
            0.95 => 0.07516774f0)
        )
    μ_200_vals = Dict(
        :ButlerVolmer=> [0.0521894f0, 0.0421092f0, 0.1289288f0], 
        :Marcus=> [0.0544722f0, 0.0466285f0, 0.2127188f0],
        :AsymptoticMarcusHushChidsey=>[0.055173985f0, 0.04740624957f0, 0.17594977f0]
        )
    μ_100_T400_vals = Dict(
        :ButlerVolmer=>[0.026958f0, 0.032575f0, 0.1537748f0] , 
        :Marcus=>[0.021007f0, 0.022126f0, 0.13133968f0],
        :AsymptoticMarcusHushChidsey=>[0.019907f0, 0.020134f0, 0.1077556f0]
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
            # TODO: fix this test for Marcus – currently, it jumps past inverted region for vector case but not for broadcast case...generally would be good to have more robust handling of cases when there are multiple valid solutions for `overpotential`...
            if !(km isa Marcus)
                @test all(isapprox.(μ_100_T400(xs), μ_100_T400.(xs), atol=5f-5)) && all(isapprox.(µ_100_T400(xs), μ_100_T400_vals[typeof(km).name.name], atol=5f-5))
            end
        end
    end
end

@testset "Kinetic g" begin
    g_200_T400_vals = Dict(
        :ButlerVolmer => [0.0211686, 0.041466, 0.079388], 
        :Marcus => [0.020072, 0.034498, 0.061939],
        :AsymptoticMarcusHushChidsey => [0.019862, 0.0331067, 0.055262]
        )
    g_50_2_vals = Dict(:ButlerVolmer=>0.03515968, :Marcus=>0.035497, :AsymptoticMarcusHushChidsey=>0.0356339)
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
    v1_vals = Dict(:ButlerVolmer=>[0.045547, 0.841501],
                :Marcus=>[0.04901, 0.81884],
                :AsymptoticMarcusHushChidsey=>[0.04958, 0.81897])
    v2_vals = Dict(:ButlerVolmer=>[0.090478, 0.77212],
                :Marcus=>[0.078934, 0.804275],
                :AsymptoticMarcusHushChidsey=>[0.07609, 0.81652])
    v3_vals = Dict(:ButlerVolmer=>[0.105747, 0.6809],
                :Marcus=>[0.154934, 0.54218],
                :AsymptoticMarcusHushChidsey=>[0.147034,0.566617])
    I_max_T300 = Dict(:ButlerVolmer=>700.0,
                :Marcus=>600.0,
                :AsymptoticMarcusHushChidsey=>650.0)
    I_max_T400 = Dict(:ButlerVolmer=>250.0,
                :Marcus=>400.0,
                :AsymptoticMarcusHushChidsey=>500.0)
    for km in kms
        model_type = nameof(typeof(km))
        @testset "$model_type" begin
            # simplest case, just one pair of x values (this function is still pretty slow though)
            v1 = find_phase_boundaries(100, km)

            @test all(isapprox.(common_tangent(v1, 100, km), Ref(0.0), atol=1e-6))
            v2 = find_phase_boundaries(100, km, T=350)     
            @test all(isapprox.(common_tangent(v2, 100, km, T=350), Ref(0.0), atol=1e-5))
            # they should get "narrower" with temperature
            @test v2[1] > v1[1]
            @test v2[2] < v1[2]
            v3 = find_phase_boundaries(400, km, guess=[0.1,0.8])
            # ...and also with current
            @test v3[1] > v1[1]
            @test v3[2] < v1[2]

            # test actual numerical values too
            @test all(isapprox.(v1, v1_vals[model_type], atol=1e-5))
            @test all(isapprox.(v2, v2_vals[model_type], atol=1e-4))
            @test all(isapprox.(v3, v3_vals[model_type], atol=1e-4))

            # and actually build the full map for a couple temps
            pbs, I = phase_diagram(km, I_max=700, I_step=50)
            @test maximum(I) == I_max_T300[model_type]
        end
    end

    @testset "MHC" begin
        # do some lighterweight stuff here since this takes awhile
        mhc = MarcusHushChidsey(amhc.A, amhc.λ)
        # also give it a good initial guess so it only takes one or two optimizer steps
        v1 = find_phase_boundaries(100, mhc, guess=v1_vals[:AsymptoticMarcusHushChidsey])
        @test all(isapprox.(v1, [0.04967204036, 0.81822937543], atol=1e-5))

        v2 = find_phase_boundaries(100, mhc, T=350, guess=v2_vals[:AsymptoticMarcusHushChidsey])
        @test all(isapprox.(v2, [0.075615, 0.8171636], atol=1e-5))

        v3 = find_phase_boundaries(400, mhc, guess=v3_vals[:AsymptoticMarcusHushChidsey])
        @test all(isapprox.(v3, [0.1467487, 0.5704125], atol=1e-5))
    end
end
