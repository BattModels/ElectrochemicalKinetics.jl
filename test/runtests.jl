using ElectrochemicalKinetics
using QuadGK
using Test

function test_vector_voltages(model, voltages)
    @test isapprox(rate_constant(voltages, model), rate_constant.(voltages, Ref(model)))
    @test isapprox(rate_constant(voltages, model, T=1000), rate_constant.(voltages, Ref(model), T=1000))
end

@testset "ElectrochemicalKinetics.jl" begin
    global Vs = 0.01:0.01:0.1
    @testset "Rate constant computation" begin
        include("non_integral_model_tests.jl")
        include("integral_model_tests.jl")
    end

    @testset "Fitting" begin
        include("fitting_tests.jl")
    end

    @testset "Phase Diagrams" begin
        include("phase_diagram_tests.jl")
    end

    @testset "Integration" begin
        include("integration.jl")
    end
end
