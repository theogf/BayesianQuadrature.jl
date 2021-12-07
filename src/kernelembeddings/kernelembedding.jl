abstract type AbstractKernelEmbedding end

measure(ke::AbstractKernelEmbedding) = ke.measure

struct KernelEmbedding{K::Kernel,Tm,Tσ::Real,Tl} <: AbstractKernelEmbedding
    kernel::K # Kernel function
    σ::Tσ # Kernel variance
    l::Tl # Kernel lengthscale
    measure::Tm # Measure
end

function KernelEmbedding(bquad::BayesQuad, prior)
    return KernelEmbedding(kernel(bquad), bquad.σ, bquad.l, prior)
end

scale(ke::KernelEmbedding) = ke.σ

@doc raw"""
    kernel_mean(ke::KernelEmbedding{Kernel,Measure}, samples::AbstractVector)

Compute the kernel mean of the kernel embedding `ke` for each one of 
the `samples` $$x_i$$:
```math
    z_i = \int k(x, x_i)d\mu(x)
```
"""
kernel_mean

function kernel_mean(ke::KernelEmbedding{<:SqExponentialKernel,<:AbstractMvNormal}, samples::AbstractVector)
    z = samples .- Ref(mean(measure(ke)))
    B = Λ(ke.l)
    return scale(ke) / sqrt(det(inv(B) * cov(measure(ke)) + I)) *
           exp.(- PDMats.invquad.(Ref(PDMat(B + cov(measure(ke)))), z) / 2)
end

@doc raw"""
    kernel_variance(ke::KernelEmbedding{Kernel,Measure})

Compute the kernel variance of the given kernel embedding:
```math
    C = \int\int k(x,x')d\mu(x)d\mu(x')
```
"""
kernel_variance


function kernel_variance(ke::KernelEmbedding{<:SqExponentialKernel,<:AbstractMvNormal})
    B = Λ(ke.l)
    return scale(ke) / sqrt(det(2 * inv(B) * cov(measure(ke)) + I))
end


# Turn the lengthscale into a Diagonal matrix of noise
Λ(l::Real) = abs2(l) * I
Λ(l::AbstractVector) = Diagonal(abs2.(l))
Λ(l::LowerTriangular) = Cholesky(l, :L, 1)
Λ(l::AbstractMatrix) = l * l'

# Turns a transform into a lengthscale
transform_param(t::ScaleTransform) = inv(first(t.s))
transform_param(t::ARDTransform) = inv.(t.v)
transform_param(t::LinearTransform) = inv(t.A)
