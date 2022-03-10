using FastGaussQuadrature: gausslegendre

# helper function to be able to "naturally" use the interval of integration that we want with FastGaussQuadrature, which returns nodes and weights for the interval [-1, 1] (see https://en.wikipedia.org/wiki/Gaussian_quadrature#Change_of_interval)
function setup_integration(f, N)
    nodes, weights = f(N)
    function scale(lb::Real, ub::Real, nodes = nodes, weights = weights)
        # use nodes and weights here
        w = @. 0.5 * (ub - lb) * weights
        x = @. 0.5 * (ub + lb) + 0.5 * (ub - lb) * nodes
        x, w
    end
    function scale(lb::T, ub::T, nodes = nodes, weights = weights) where T <: AbstractVector
        ns_and_ws = scale.(lb, ub, Ref(nodes), Ref(weights))
        ns, ws = Zygote.unzip(ns_and_ws)::Vector{Tuple{T, T}}
        reduce(hcat, ns), reduce(hcat, ws)
    end
end

const scale = setup_integration(gausslegendre, 1000)

# const nodes, weights = gausslegendre(1000)
# function scale(lb::Real, ub::Real, nodes = nodes, weights = weights)
#     # use nodes and weights here
#     w = @. 0.5 * (ub - lb) * weights
#     x = @. 0.5 * (ub + lb) + 0.5 * (ub - lb) * nodes
#     x, w
# end
# function scale(lb::AbstractVector, ub::AbstractVector, nodes::AbstractVector = nodes, weights::AbstractVector = weights)
#     ns_and_ws = scale.(lb, ub, Ref(nodes), Ref(weights))
#     @show typeof(ns_and_ws)
#     ns, ws = Zygote.unzip(ns_and_ws)::Vector{Tuple{Vector{Float64}, Vector{Float64}}}
#     reduce(hcat, ns), reduce(hcat, ws)
# end

# Can't get this to broadcast correctly - but this is the right sort of
# approach to take here. Of note - this is 2x faster in the forwards pass
# and about 10x faster in the backwards pass compared to the other method
# function scale_integration_nodes2(unscaled_x, unscaled_w, lb, ub)
#     w_scale = @. 0.5 * (ub - lb)
#     x_scale = @. 0.5 * (ub + lb) + 0.5 * (ub - lb)
#     unscaled_x * x_scale', unscaled_w * w_scale'
#     # x_scale * unscaled_x, unscaled_w * w_scale'
# end
