using ElectrochemicalKinetics
using Test

@testset "ElectrochemicalKinetics.jl" begin
    @testset "Rate constant computation" begin
        include("non_integral_model_tests.jl")
        include("integral_model_tests.jl")
    end

    @testset "Fitting" begin
        include("fitting_tests.jl")
    end
end
