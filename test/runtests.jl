using ElectrochemicalKinetics
using Test

@testset "ElectrochemicalKinetics.jl" begin
    global Vs = 0.01:0.01:0.1
    @testset "Rate constant computation" begin
        include("non_integral_model_tests.jl")
        include("integral_model_tests.jl")
    end

    @testset "Fitting" begin
        include("fitting_tests.jl")
    end

    @testset "Thermo" begin
        include("thermo_tests.jl")
    end
end
