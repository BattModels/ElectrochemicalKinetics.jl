# test values come either from analytical inversion of expression (BV) or references evaluated as of 2/22/22 and presumed to be correct

# TODO: throughout, add tests of varying kT

@testset "Non-Integral Models" begin
    @testset "Butler-Volmer" begin
        bv1 = ButlerVolmer()
        bv2 = ButlerVolmer(3.0)
        bv3 = ButlerVolmer(1.0, 0.2)

        @testset "Scalars" begin
            # zero voltage
            @test all(rate_constant.(Ref(0), [bv1, bv2, bv3]).==0)

            # oxidative
            for bv in [bv1, bv2, bv3]
                @test bv(.026/bv.α, true) ≈ bv.A * ℯ
            end

            # reductive
            for bv in [bv1, bv2, bv3]
                @test bv(-.026/(1-bv.α), false) ≈ bv.A * ℯ
            end

            # net
            for bv in [bv1, bv2, bv3]
                @test bv(.026/bv.α) ≈ bv.A * (ℯ - exp((bv.α-1)/bv.α))
            end
        end

        @testset "Vectors" begin
            @testset "Vector Voltages" begin
                test_vector_voltages(bv1, Vs)
                test_vector_voltages(bv2, Vs)
                test_vector_voltages(bv3, Vs)
                for model in [bv1, bv2] # bv3 won't be symmetrical
                    @test isapprox(model(Vs), model(-Vs)) #symmetry
                end
            end
            @testset "Vector Models" begin
                A_vals = [1.0, 2.0]
                α_vals = [0.4, 0.5, 0.6]
                bvv1 = ButlerVolmer(A_vals) # autopopulate α, combine scalar/vector
                bvv2 = ButlerVolmer(A_vals, α_vals[1:2]) # both vectors
                @test_throws DimensionMismatch ButlerVolmer(A_vals, α_vals) # invalid, different lengths

                # test that it's equivalent to the scalar model version
                @test rate_constant(0.1, bvv1) == rate_constant.(Ref(0.1), ButlerVolmer.(A_vals))
                @test rate_constant(0.1, bvv2) == rate_constant.(Ref(0.1), ButlerVolmer.(A_vals, α_vals[1:2]))
            end
            @testset "Vector Voltages and Models" begin
                # TODO: this
            end
        end
    end

    @testset "Marcus" begin
        m1 = Marcus(0.25)
        m2 = Marcus(10, 0.25)
        @testset "Scalars" begin
            # TODO: this could all be more thorough, I literally just picked some arbitrary points
            @test m1(0.0) == 0.0
            @test isapprox(m1(0.05, true), 0.031381, atol=1e-6)
            @test isapprox(m1(0.2, false), 0.908324, atol=1e-6)
            @test isapprox(m1(0.13), 0.570862, atol=1e-6)
            @test isapprox(m2(0.13), 5.708625, atol=1e-6)
        end

        @testset "Vectors" begin
            @testset "Vector Voltages" begin
                test_vector_voltages(m1, Vs)
                test_vector_voltages(m2, Vs)
                for model in [m1, m2]
                    @test isapprox(model(Vs), model(-Vs)) #symmetry
                end
            end
            @testset "Vector Models" begin
                A_vals = [1, 100, 10000]
                λ_vals = [0.1, 0.3, 0.7]
                mv1 = Marcus(λ_vals)
                mv2 = Marcus(A_vals, λ_vals)

                @test rate_constant(0.1, mv1) == rate_constant.(Ref(0.1), Marcus.(λ_vals))
                @test rate_constant(0.1, mv2) == rate_constant.(Ref(0.1), Marcus.(A_vals, λ_vals))
            end
            @testset "Vector Voltages and Models" begin
                # TODO: this
            end
        end
    end

    @testset "asymptotic MHC" begin
        amhc = AsymptoticMarcusHushChidsey(0.25)
        @testset "Scalars" begin
            @test amhc(0.0) == 0.0
            @test amhc(0.0, true) == amhc(0.0, false) ≈ 0.00596440064
            @test amhc(0.2) == amhc(-0.2) ≈ 0.1006330238
            @test amhc(0.2, kT=.015) ≈ 0.06360960459
            @test amhc(0.2, kT=.03) ≈ 0.11255129857   
        end
        @testset "Vector Voltages" begin
            test_vector_voltages(amhc, Vs)
            @test isapprox(amhc(Vs), amhc(-Vs)) #symmetry
        end
        @testset "Vector Models" begin
            A_vals = [1, 100, 10000]
            λ_vals = [0.1, 0.3, 0.7]
            amhc1 = AsymptoticMarcusHushChidsey(λ_vals)
            amhc2 = AsymptoticMarcusHushChidsey(A_vals, λ_vals)

            @test rate_constant(0.1, amhc1) == rate_constant.(Ref(0.1), AsymptoticMarcusHushChidsey.(λ_vals))
            @test rate_constant(0.1, amhc2) == rate_constant.(Ref(0.1), AsymptoticMarcusHushChidsey.(A_vals, λ_vals))
        end
        @testset "Vector Voltages and Models" begin
            # TODO: this
        end
    end
end
