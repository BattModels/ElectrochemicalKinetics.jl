using ElectrochemicalKinetics
using QuadGK
using Test
using Unitful
import Unitful: @u_str

function test_vector_voltages(model, voltages)
    @test isapprox(rate_constant(voltages, model), rate_constant.(voltages, Ref(model)))
    @test isapprox(rate_constant(voltages, model, T=1000), rate_constant.(voltages, Ref(model), T=1000))
end

function test_vector_models(ModelType, params)
    model1 = ModelType(params[2:end]...)
    model2 = ModelType(params...)
    broadcast_params = []
    for p in params
        if length(p) == 1
            push!(broadcast_params, Ref(p))
        else
            push!(broadcast_params, p)
        end
    end
    @test isapprox(rate_constant(0.1, model1), rate_constant.(Ref(0.1), ModelType.(broadcast_params[2:end]...)))
    @test isapprox(rate_constant(0.1, model1, T=1000), rate_constant.(Ref(0.1), ModelType.(broadcast_params[2:end]...), T=1000))
    @test isapprox(rate_constant(0.1, model2), rate_constant.(Ref(0.1), ModelType.(broadcast_params...)))
    @test isapprox(rate_constant(0.1, model2, T=200), rate_constant.(Ref(0.1), ModelType.(broadcast_params...), T=200))
    #with units
    @test isapprox(rate_constant(0.1u"V", model1), rate_constant.(Ref(0.1), ModelType.(broadcast_params[2:end]...)))
    @test isapprox(rate_constant(0.1u"V", model1, T=1000u"K"), rate_constant.(Ref(0.1), ModelType.(broadcast_params[2:end]...), T=1000))
    @test isapprox(rate_constant(0.1u"V", model2), rate_constant.(Ref(0.1), ModelType.(broadcast_params...)))
    @test isapprox(rate_constant(0.1u"V", model2, T=200u"K"), rate_constant.(Ref(0.1), ModelType.(broadcast_params...), T=200))

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
