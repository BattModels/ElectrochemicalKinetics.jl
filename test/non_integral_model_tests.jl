# make sure to test each one for:
# - single V value, ox, red, and net
# - vector of V values, ox, red, and net

@testset "Non-Integral Models" begin
    @testset "Butler-Volmer" begin
        bv1 = ButlerVolmer()
    end

    @testset "Marcus" begin
        
    end

    @testset "asymptotic MHC" begin
        
    end
end
