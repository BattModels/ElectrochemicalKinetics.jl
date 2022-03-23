@testset "Integration" begin
  lb, ub = 0.1, 0.2
  nodes, weights = ElectrochemicalKinetics.scale(lb, ub)
  @test nodes[1] ≈ lb atol = 1e-5
  @test nodes[end] ≈ ub atol = 1e-5

  lb = rand(2)
  ub = lb .+ rand.()

  nodes, weights = ElectrochemicalKinetics.scale(lb, ub)
  @test nodes isa Matrix
  @test nodes[1, :] ≈ lb atol = 1e-5
  @test nodes[end, :] ≈ ub atol = 1e-5
end
