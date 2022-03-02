# test values come either from analytical inversion of expression (BV) or references evaluated as of 2/22/22 and presumed to be correct

@testset "Non-Integral Models" begin
    @testset "Butler-Volmer" begin
        bv1 = ButlerVolmer()
        bv2 = ButlerVolmer(3.0)
        bv3 = ButlerVolmer(1.0, 0.2)

        @testset "Scalar Arguments" begin
            # zero voltage
            @test all(compute_k.(Ref(0), [bv1, bv2, bv3]).==0)

            # oxidative
            for bv in [bv1, bv2, bv3]
                @test bv(.026/bv.α, Val(true)) ≈ bv.A * ℯ
            end

            # reductive
            for bv in [bv1, bv2, bv3]
                @test bv(-.026/(1-bv.α), Val(false)) ≈ bv.A * ℯ
            end

            # net
            for bv in [bv1, bv2, bv3]
                @test bv(.026/bv.α) ≈ bv.A * (ℯ - exp((bv.α-1)/bv.α))
            end
        end

        @testset "Vector Arguments" begin
            @test bv1(.026 .* collect(0:2:6), Val(true)) == exp.(0:3)
            @test bv1(.026 .* collect(0:2:6)) == exp.(0:3) - exp.(-(0:3))
            # TODO: test bv2, bv3
        end
    end

    @testset "Marcus" begin
        # TODO: this could all be more thorough, I literally just picked some arbitrary points
        m = Marcus(0.25)
        @test m(0.0) == 0.0
        @test isapprox(m(0.05, Val(true)), 0.031381, atol=1e-6)
        @test isapprox(m(0.2, Val(false)), 0.908324, atol=1e-6)
        @test isapprox(m(0.13), 0.570862, atol=1e-6)
        @test m(0.13) == m(-0.13)

        # test vector voltages and also symmetry
        @test m(Vs) == m(-Vs)
    end

    @testset "asymptotic MHC" begin
        amhc = AsymptoticMarcusHushChidsey(0.25)
        @test amhc(0.0) == 0.0
        @test amhc(0.0, true) == amhc(0.0, false) ≈ 0.00596440064
        @test amhc(0.2) == amhc(-0.2) ≈ 0.1006330238
        @test amhc(0.2, kT=.015) ≈ 0.06360960459
        @test amhc(0.2, kT=.03) ≈ 0.11255129857

        vec_out = amhc(Vs)
        @test length(vec_out) == length(Vs)
        @test vec_out[1] == amhc(Vs[1])
    end
end
