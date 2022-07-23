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
        ns, ws = Zygote.unzip(ns_and_ws)::Tuple{Vector{T}, Vector{T}}
        reduce(hcat, ns), reduce(hcat, ws)
    end
end

const scale = setup_integration(gausslegendre, 500)
const scale_coarse = setup_integration(gausslegendre, 50)
