"""
    BayesQuad(k::Kernel; l=1.0, σ::Real=1.0)

Bayesian Quadrature object.
You can pass any kernel and the lengthscale and variance will be extracted.
`l` can be a `Real`, a `AbstractVector` or a `LowerTriangular`.
"""
struct BayesQuad{TK,Tl,Tσ} <: AbstractBQ
    kernel::TK
    l::Tl
    σ::Tσ
end

function BayesQuad(k::Kernel; l=1.0, σ::Real=1.0)
    σ > 0 || ArgumentError("σ should be positive")
    l isa AbstractMatrix && (
        l isa LowerTriangular || throw(
            ArgumentError(
                "For l an AbstractMatrix, only LowerTriangular matrices are accepted"
            ),
        )
    )
    return BayesQuad(k, l, σ)
end

function BayesQuad(k::TransformedKernel{TK,Tt}; l=nothing, σ=1.0) where {TK,Tt}
    Tt <: Union{ScaleTransform,ARDTransform,LinearTransform} ||
        error("No lengthscale could be extracted from kernel $k,\n
        only ScaleTransform, ARDTransform and LinearTransform are allowed")
    l = param(k.transform)

    return BayesQuad(k.kernel; l=l, σ=σ)
end

function BayesQuad(k::ScaledKernel; l=1.0, σ=nothing)
    σ = first(k.σ²)
    return BayesQuad(k.kernel; l=l, σ=σ)
end

function kernel(b::BayesQuad)
    return b.σ * transform(b.kernel, inv.(b.l))
end

function kernel(b::BayesQuad{<:Kernel,<:LowerTriangular})
    return b.σ * transform(b.kernel, LinearTransform(inv(b.l)))
end

function quadrature(
    bquad::BayesQuad{<:SqExponentialKernel},
    model::AbstractBQModel{<:MvNormal},
    samples,
)
    isempty(samples) && error("The collection of samples is empty")
    y = integrand(model).(samples)
    K = kernelpdmat(kernel(bquad), samples)
    z = calc_z(samples, prior(model), bquad)
    C = calc_C(prior(model), bquad)
    var = evaluate_var(z, K, C)
    return Normal(evaluate_mean(z, K, y), max(var,zero(var)))
end

Λ(bquad::BayesQuad{<:SqExponentialKernel,<:Real}) = abs2(bquad.l) * I
Λ(bquad::BayesQuad{<:SqExponentialKernel,<:AbstractVector}) = Diagonal(abs2.(bquad.l))
Λ(bquad::BayesQuad{<:SqExponentialKernel,<:LowerTriangular}) = bquad.l * bquad.l'

scale(bquad::BayesQuad) = bquad.σ

function evaluate_mean(z, K, y)
    return dot(z, K \ y)
end

function evaluate_var(z, K, C)
    return C - PDMats.invquad(K, z)
end
