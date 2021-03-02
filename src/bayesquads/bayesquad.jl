"""
    BayesQuad(k::Kernel; l=1.0, σ::Real=1.0)

Bayesian Quadrature object.
You can pass any kernel and the lengthscale and variance will be extracted.
`l` can be a `Real`, a `AbstractVector` or a `LowerTriangular`.
"""
struct BayesQuad{TK,Tl,Tσ} <: AbstractBayesQuad{TK, Tl}
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

function quadrature(
    bquad::BayesQuad{<:SqExponentialKernel},
    model::AbstractBayesQuadModel{<:MvNormal},
    samples,
)
    isempty(samples) && error("The collection of samples is empty")
    y = integrand(model).(samples)
    K = kernelpdmat(kernel(bquad), samples)
    z = calc_z(samples, prior(model), bquad)
    C = calc_C(prior(model), bquad)
    return Normal(evaluate_mean(z, K, y), sqrt(evaluate_var(z, K, C)))
end
