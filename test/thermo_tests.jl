include("../examples/thermo.jl")

@testset "Basic free energy functions" begin
    # enthalpy of mixing
    @test hs(0.5) == 0.025 # value
    @test hs(0.1) ≈ hs(0.9) # symmetry
    @test hs([0.1, 0.2, 0.7]) ≈ [0.009, 0.016, 0.021] # vector input
    @test hs(0.5; Ω=1) == 0.25 # kwargs

    # similar set for entropy...
    @test s(0.5) ≈ -kB * 2 * 0.5 * log(0.5)
    @test s(0.1) == s(0.9)
    @test s([0.1, 0.2, 0.7]) ≈ [2.8012399817e-5, 4.3119676836e-5, 5.2638176908e-5]

    # ...and thermodynamic free energy...
    @test g_thermo(0.5) == hs(0.5) - T * s(0.5) + muoA * 0.5 + muoB * 0.5 ≈ 0.032200909
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

xs = [0.1, 0.5, 0.95]
km = ButlerVolmer()
μ_50 = μ_kinetic(50, km)
μ_200 = μ_kinetic(200, km)
μ_100_T400 = μ_kinetic(100, km, T=400)

@testset "Kinetic μ" begin
    # test that the right (approximate) relationships hold
    μs = Dict(0.1 => 0.24249555384, 0.5 => 0.249469635, 0.95 => 0.35480771)
    for x in xs
        @test μ_50(x) ≈ μs[x] ≈ μ_thermo(x) + fit_overpotential((1-x)*km, 50)
    end

    # test vector inputs
    @test_throws MethodError μ_200(xs)
    @test μ_200.(xs) ≈ [0.314567655, 0.32155307, 0.4268963]
    @test μ_100_T400.(xs) ≈ [0.259213284, 0.285511025, 0.41673288]
end

@testset "Kinetic g" begin
    g_50 = g_kinetic(50, km)
    g_200_T400 = g_kinetic(200, km, T=400)

    # integral of the derivative should be the original fcn (up to a constant, which we know)
    integrated_μs = [v[1] for v in quadgk.(μ_50, zero(xs), xs)]
    @test g_50(xs) ≈ integrated_μs .+ muoA

    # check scalar input works
    @test g_50(xs[2]) == g_50(xs)[2]

    # check a few other values, these are just from  me trusting that the function is running correctly when I write these...
    @test g_200_T400(xs) ≈ [0.04661525, 0.17184189, 0.33075273]
end
