"""
    BayesQuad(k::Kernel; l=1.0, σ::Real=1.0)

Tool for estimating the bayesian quadrature for a given `BayesQuadModel` and a `Sampler`.

`BayesQuad` estimate the probability distribution `p(I)` of I = ∫ f(x) p(x) dx
## Argument
- `k::Kernel` : kernel from `KernelFunctions.jl`, the kernel cannot be composite
## Keywords argument
- `l` : lengthscale of the kernel. It can be:
    - `Real` : isotropic kernel
    - `AbstractVector` : ARD kernel (one lengthscale per dimension)
    - `LowerTriangular` : Linear transformation of the inputs
- `σ` : variance of the kernel (`k(x,x') -> σ * k(x,x')`)

Note: If `k` already has a variance and/or a transformation,
these will be automatically extracted
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
