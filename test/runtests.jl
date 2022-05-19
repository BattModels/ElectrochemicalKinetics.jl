using ElectrochemicalKinetics
using QuadGK
using Test

function test_vector_voltages(model, voltages)
    @test isapprox(compute_k(voltages, model), compute_k.(voltages, Ref(model)))
    @test isapprox(compute_k(voltages, model, kT=0.4), compute_k.(voltages, Ref(model), kT=0.4))
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
