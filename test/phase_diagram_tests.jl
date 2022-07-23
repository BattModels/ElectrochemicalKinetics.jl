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
        :ButlerVolmer=>[0.023721278f0, 0.02681895f0, 0.120048858f0] , 
        :Marcus=>[0.02481338f0, 0.0288534856f0, 0.1539812383f0],
        :AsymptoticMarcusHushChidsey=>[0.0252367557f0, 0.029513253f0, 0.146351765f0]
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
            @test μ_100_T400(xs) == μ_100_T400.(xs) && all(isapprox.(µ_100_T400(xs), μ_100_T400_vals[typeof(km).name.name], atol=1f-5))
        end
    end
end

@testset "Kinetic g" begin
    g_200_T400_vals = Dict(
        :ButlerVolmer => [0.020563, 0.03755, 0.0661336], 
        :Marcus => [0.0207786, 0.0390248, 0.074701],
        :AsymptoticMarcusHushChidsey => [0.0208467, 0.03939946, 0.0740394]
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
    @test all(isapprox.(v1, [0.045547, 0.841501], atol=1e-5))
    @test all(isapprox.(v2, [0.083289, 0.796802], atol=1e-4))
    @test all(isapprox.(v3, [0.105747, 0.6809], atol=1e-4))
end
