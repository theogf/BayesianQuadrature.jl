include("bayesquad.jl")
include("logbayesquad.jl")

function kernel(b::AbstractBayesQuad)
    return b.σ * transform(b.kernel, inv.(b.l))
end

function kernel(b::AbstractBayesQuad{<:Kernel,<:LowerTriangular})
    return b.σ * transform(b.kernel, LinearTransform(inv(b.l)))
end

Λ(bquad::AbstractBayesQuad{<:SqExponentialKernel,<:Real}) = abs2(bquad.l) * I
Λ(bquad::AbstractBayesQuad{<:SqExponentialKernel,<:AbstractVector}) = Diagonal(abs2.(bquad.l))
Λ(bquad::AbstractBayesQuad{<:SqExponentialKernel,<:LowerTriangular}) = bquad.l * bquad.l'

scale(bquad::AbstractBayesQuad) = bquad.σ
lengthscale(bquad::AbstractBayesQuad) = bquad.l

function evaluate_mean(z, K, y)
    return dot(z, K \ y)
end

function evaluate_var(z, K, C)
    return C - PDMats.invquad(K, z)
end
