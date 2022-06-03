# all these numbers are just references evaluated as of 2/22/22
# TODO: add tests of `integrand` function
# TODO: throughout, add tests of varying kT
@testset "Integral Models" begin
    @testset "MarcusHushChidsey" begin
        mhc = MarcusHushChidsey(0.25)
        @testset "Scalars" begin
            # test some values
            @test rate_constant(0, mhc) == 0.0
            @test rate_constant(0, mhc, true) == rate_constant(0, mhc, false) ≈ 0.0061371322
            @test isapprox(rate_constant(0.1, mhc), 0.031239978, atol=1e-6)
            
            # test net rates and symmetry
            @test isapprox(rate_constant(0.1, mhc, true) - rate_constant(0.1, mhc, false), rate_constant(0.1, mhc), atol=1e-6)
            @test rate_constant(-0.1, mhc) == rate_constant(0.1, mhc)
        end

        @testset "Vector Voltages" begin
            test_vector_voltages(mhc, Vs)
            @test isapprox(mhc(Vs), mhc(-Vs)) #symmetry
            # test that it's close to asymptotic version...arguably this should perhaps be in the AMHC tests
            amhc = AsymptoticMarcusHushChidsey(0.25)
            @test all(isapprox.(rate_constant(Vs, mhc), rate_constant(Vs, amhc), atol=2e-4))
        end
        
        @testset "Vector Models" begin
            A_vals = [1, 10, 1000]
            λ_vals = [0.1, 0.4, 0.7]
            avg_dos_val = 1.0
            test_vector_models(MarcusHushChidsey, [A_vals, λ_vals, avg_dos_val])
        end
    end

    @testset "MarcusHushChidseyDOS" begin
        # test that a uniform DOS version matches MHC
        flat_dos = [-5 1; 5 1]
        mhcd = MarcusHushChidseyDOS(0.25, flat_dos)
        @test isapprox(rate_constant(0, mhcd), 0.0, atol=1e-6)
        mhc = MarcusHushChidsey(0.25)
        mhcd_k = rate_constant(Vs, mhcd)
        mhc_k = rate_constant(Vs, mhc)
        @test all(isapprox.(mhc_k, mhcd_k, atol=1e-6))

        # make some cartoon smooth DOSes...
        Es = -5:0.1:5
        gaussian(x) = exp(-x^2)
        # here's one with two overlapping peaks
        dos = hcat(Es, gaussian.(0.7 .*(Es .- 2)) .+ gaussian.(0.7 .*(Es .+ 2)))
        mhcd = MarcusHushChidseyDOS(0.25, dos)
        # this should be symmetric
        @test all(isapprox.(rate_constant(Vs, mhcd), rate_constant(-Vs, mhcd), atol=1e-6))
        # now try an asymmetric one
        dos = hcat(Es .+ 1, dos[:,2])
        mhcd = MarcusHushChidseyDOS(0.25, dos)
        @test all(rate_constant(-Vs, mhcd) .> rate_constant(Vs, mhcd))
    end

    @testset "Quantum Capacitance" begin
        # TODO: this
    end
end