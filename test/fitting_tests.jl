# TODO: add tests for negative overpotentials

@testset "Overpotential Fitting" begin
    params = Dict(ButlerVolmer=>[], Marcus=>[0.25], AsymptoticMarcusHushChidsey=>[0.25], MarcusHushChidsey=>[0.25])
    for model_type in [ButlerVolmer] #, Marcus, AsymptoticMarcusHushChidsey, MarcusHushChidsey]
        @testset "$model_type" begin
            m =  model_type(params[model_type]...)
            # test a number less than the initial guess
            V1 = 0.05
            k1 = rate_constant(V1, m)
            @test isapprox(overpotential(k1, m)[1], V1, atol=1e-5)
            # and one greater
            V2 = 0.15
            k2 = rate_constant(V2, m)
            @test isapprox(overpotential(k2, m)[1], V2, atol=1e-5)
            # test that vectorized fitting works...first for one model and multiple k's
            # this test is turned off for now because better adjoint of overpotential currently doesn't work with it...
            # @test isapprox(overpotential(m, [k1, k2]), [V1, V2], atol=5e-4)
            # TODO: add tests for vector model case
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
            k1 = rate_constant(V1, m)
            @test isapprox(overpotential(k1, m)[1], V1, atol=1e-5)
            V2 = 0.15
            k2 = rate_constant(V2, m)
            @test isapprox(overpotential(k2, m)[1], V2, atol=1e-5)
        end
    end
end

@testset "Model Parameter Fitting" begin
    
end