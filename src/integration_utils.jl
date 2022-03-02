using FastGaussQuadrature

# helper function to be able to "naturally" use the interval of integration that we want with FastGaussQuadrature, which returns nodes and weights for the interval [-1, 1] (see https://en.wikipedia.org/wiki/Gaussian_quadrature#Change_of_interval)
function get_nodes_weights(lb, ub; N=1000,  quadfun=gausslegendre)
    unscaled_x, unscaled_w = quadfun(N)
    w = 0.5 * (ub - lb) * unscaled_w
    x = 0.5 * (ub + lb) + 0.5 * (ub - lb) * unscaled_x
    w, x
end