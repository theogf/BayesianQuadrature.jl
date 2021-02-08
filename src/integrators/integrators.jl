include("bquad.jl")

gp(i::AbstractIntegrator) = i.gp

kernel(i::AbstractIntegrator) = kernel(gp(i))
kernel(gp::GP) = gp.kernel

function evaluate_mean(z, K, y)
    return dot(z, K \ y)
end

function evaluate_var(z, K, C)
    return C - invquad(K, z)
end