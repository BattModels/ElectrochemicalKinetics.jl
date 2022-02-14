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
