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
these will be automatically extracted and replace the given keyword arguments
"""
struct BayesQuad{TK,Tl,Tσ} <: AbstractBQ{TK}
    kernel::TK
    l::Tl
    σ::Tσ
end

function BayesQuad(k::Kernel; l=1.0, σ::Real=1.0)
    k, (l, σ) = get_kernel_params(k; l, σ)
    check_kernel_parameters(l, σ)
    return BayesQuad(k, l, σ)
end



function quadrature(
    bquad::BayesQuad{<:SqExponentialKernel},
    model::AbstractBQModel{<:MvNormal},
    samples,
)
    isempty(samples) && error("The collection of samples is empty")
    y = integrand(model).(samples)
    K = kernelpdmat(kernel(bquad), samples)
    z = calc_z(samples, p_0(model), bquad)
    C = calc_C(p_0(model), bquad)
    var = evaluate_var(z, K, C)
    if var < 0
        if var > -1e-5
            @warn "Variance was negative (numerical error) and set to 0"
        else
            error("Obtained variance was negative")
        end
    end
    return Normal(evaluate_mean(z, K, y), max(var, zero(var)))
end

scale(bquad::BayesQuad) = bquad.σ

function evaluate_mean(z, K, y)
    return dot(z, K \ y)
end

function evaluate_var(z, K, C)
    return C - PDMats.invquad(K, z)
end
