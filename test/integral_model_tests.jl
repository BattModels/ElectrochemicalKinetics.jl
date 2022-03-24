# all these numbers are just references evaluated as of 2/22/22

@testset "Integral Models" begin
    @testset "MHC" begin
        mhc = MarcusHushChidsey(0.25)

        # test some values
        @test compute_k(0, mhc) == 0.0
        @test compute_k(0, mhc, true) == compute_k(0, mhc, false) ≈ 0.0061371322
        @test isapprox(compute_k(0.1, mhc), 0.031239978, atol=1e-6)
        
        # test net rates and symmetry
        @test isapprox(compute_k(0.1, mhc, true) - compute_k(0.1, mhc, false), compute_k(0.1, mhc), atol=1e-6)
        @test compute_k(-0.1, mhc) == compute_k(0.1, mhc)

        # test that it's close to asymptotic version, also that vector inputs work appropriately
        amhc = AsymptoticMarcusHushChidsey(0.25)
        @test all(isapprox.(compute_k(Vs, mhc), compute_k(Vs, amhc), atol=2e-4))
    end

    @testset "MHC+DOS" begin
        # test that a uniform DOS version matches MHC
        flat_dos = [-5 1; 5 1]
        mhcd = MarcusHushChidseyDOS(0.25, flat_dos)
        @test isapprox(compute_k(0, mhcd), 0.0, atol=1e-6)
        mhc = MarcusHushChidsey(0.25)
        mhcd_k = compute_k(Vs, mhcd)
        mhc_k = compute_k(Vs, mhc)
        @test all(isapprox.(mhc_k, mhcd_k, atol=1e-6))

        # make some cartoon smooth DOSes...
        Es = -5:0.1:5
        gaussian(x) = exp(-x^2)
        # here's one with two overlapping peaks
        dos = hcat(Es, gaussian.(0.7 .*(Es .- 2)) .+ gaussian.(0.7 .*(Es .+ 2)))
        mhcd = MarcusHushChidseyDOS(0.25, dos)
        # this should be symmetric
        @test all(isapprox.(compute_k(Vs, mhcd), compute_k(-Vs, mhcd), atol=1e-6))
        # now try an asymmetric one
        dos = hcat(Es .+ 1, dos[:,2])
        mhcd = MarcusHushChidseyDOS(0.25, dos)
        @test all(compute_k(-Vs, mhcd) .> compute_k(Vs, mhcd))
    end

    @testset "Quantum Capacitance" begin
        # TODO: this
    end
end