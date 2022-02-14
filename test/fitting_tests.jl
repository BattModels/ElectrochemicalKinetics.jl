@testset "Overpotential Fitting" begin
    params = Dict(ButlerVolmer=>[], Marcus=>[0.25], AsymptoticMarcusHushChidsey=>[0.25], MarcusHushChidsey=>[0.25])
    for model_type in [ButlerVolmer, Marcus, AsymptoticMarcusHushChidsey, MarcusHushChidsey]
        @testset "$model_type" begin
            m =  model_type(params[model_type]...)
            # test a number less than the initial guess
            V1 = 0.05
            k1 = compute_k(V1, m)
            @test isapprox(fit_overpotential(m, k1)[1], V1, atol=1e-5)
            # and one greater
            V2 = 0.15
            k2 = compute_k(V2, m)
            @test isapprox(fit_overpotential(m, k2)[1], V2, atol=1e-5)
            # test that vectorized fitting works...first for one model and multiple k's
            # (a bit weird that this doesn't give exactly equal results to the single-variable case, but it does get reasonably close...)
            @test isapprox(fit_overpotential(m, [k1, k2]), [V1, V2], atol=5e-4)
            # now test for one k and multiple models
            @test isapprox(fit_overpotential([m, m], k1), [V1, V1], atol=5e-4)
            # make sure that both fails since indexing would be ambiguous
            @test_throws MethodError fit_overpotential([m, m], [k1, k2])
        end
    end
    # for now, writing these separately to test different input DOSes
    @testset "MarcusHushChidseyDOS" begin
        flat_dos = [-5 1; 5 1]
        Es = -5:0.1:5
        gaussian(x) = exp(-x^2)
        # here's one with two overlapping peaks
        smooth_dos = hcat(Es, gaussian.(0.7 .*(Es .- 2)) .+ gaussian.(0.7 .*(Es .+ 2)))
        # real_dos = ...
        # for dos in [flat_dos, smooth_dos, real_dos]
        for dos in [flat_dos, smooth_dos]
            m = MarcusHushChidseyDOS(0.25, dos)
            V1 = 0.05
            k1 = compute_k(V1, m)
            @test isapprox(fit_overpotential(m, k1)[1], V1, atol=1e-5)
            V2 = 0.15
            k2 = compute_k(V2, m)
            @test isapprox(fit_overpotential(m, k2)[1], V2, atol=1e-5)
            @test isapprox(fit_overpotential(m, [k1, k2]), [V1, V2], atol=5e-4)
        end
    end
end

@testset "Model Parameter Fitting" begin
    
end