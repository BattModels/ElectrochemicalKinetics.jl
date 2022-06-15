using ElectrochemicalKinetics
using QuadGK
using Test

function test_vector_voltages(model, voltages)
    @test isapprox(rate_constant(voltages, model), rate_constant.(voltages, Ref(model)))
    @test isapprox(rate_constant(voltages, model, kT=0.1), rate_constant.(voltages, Ref(model), kT=0.1))
end

function test_vector_models(ModelType, params)
    model1 = ModelType(params[2:end]...)
    model2 = ModelType(params...)
    # @test isapprox(rate_constant(0.1, model1), rate_constant.(Ref(0.1), ModelType.(ElectrochemicalKinetics.iterate_props(model1, exclude_props = [:A]))))
    broadcast_params = []
    for p in params
        if length(p) == 1
            push!(broadcast_params, Ref(p))
        else
            push!(broadcast_params, p)
        end
    end
    @test isapprox(rate_constant(0.1, model1), rate_constant.(Ref(0.1), ModelType.(broadcast_params[2:end]...)))
    @test isapprox(rate_constant(0.1, model1, kT=0.1), rate_constant.(Ref(0.1), ModelType.(broadcast_params[2:end]...), kT=0.1))
    @test isapprox(rate_constant(0.1, model2), rate_constant.(Ref(0.1), ModelType.(broadcast_params...)))
    @test isapprox(rate_constant(0.1, model2, kT=0.1), rate_constant.(Ref(0.1), ModelType.(broadcast_params...), kT=0.1))
end

function test_vector_both(ModelType, params, voltages)
    vec_model = ModelType(params...)
    model_vec = [ModelType(p...) for p in ElectrochemicalKinetics.iterate_props(vec_model)]
    res1 = rate_constant(voltages, vec_model)
    res2 = map(t->rate_constant(t...), Iterators.product(voltages, model_vec))
    @test all(res1 .== res2)
    @test size(res1) == (length(voltages), length(vec_model))
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
