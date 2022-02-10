@testset begin "Overpotential Fitting" begin
    bv =  ButlerVolmer()
    # test a number less than the initial guess
    V1 = 0.05
    k1 = bv(V1)
    @test isapprox(fit_overpotential(bv, k1)[1], V1, atol=1e-5)
    # and one greater
    V2 = 0.15
    k2 = bv(V2)
    @test isapprox(fit_overpotential(bv, k2)[1], V2, atol=1e-5)
    # test that vectorized fitting works...this test currently fails
    @test isapprox(fit_overpotential(bv, [k1, k2]), [V1, V2], atol=1e-5)
end

@testset begin "Model Parameter Fitting" begin
    
end