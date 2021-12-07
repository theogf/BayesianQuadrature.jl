include("bayesquad.jl")
include("logbayesquad.jl")

function check_kernel_parameters(l, σ)
    σ > 0 || throw(ArgumentError("σ should be positive"))
    l isa AbstractMatrix && (
        l isa LowerTriangular || throw(
            ArgumentError(
                "For l an AbstractMatrix, only LowerTriangular matrices are accepted"
            ),
        )
    )
    l isa Union{AbstractVector,Real} && (
        all(>=(0), l) || throw(ArgumentError("At least one lengthscale value is negative"))
    )
    return nothing
end

function get_kernel_params(k::Kernel; kwargs...)
    return (k, (; kwargs...))
end

function get_kernel_params(k::ScaledKernel; kwargs...)
    # return get_kernel_params(k.kernel; merge((;kwargs...), (σ=first(k.σ²),)))
    return get_kernel_params(k.kernel; kwargs..., σ=first(k.σ²))
end

function get_kernel_params(k::TransformedKernel; kwargs...)
    check_transform(k.transform)
    return get_kernel_params(k.kernel; kwargs..., l=transform_param(k.transform))
end

function check_transform(transform)
    transform isa Union{ScaleTransform,ARDTransform,LinearTransform} ||
        throw(ArgumentError("No lengthscale could be extracted from kernel,\n
        only ScaleTransform, ARDTransform and LinearTransform are allowed"))
    return nothing
end

function kernel(b::AbstractBQ)
    if b.l isa AbstractMatrix
        b.σ * (b.kernel ∘ LinearTransform(inv(b.l)))
    else
        return b.σ * (b.kernel ∘ ScaleTransform(inv.(b.l)))
    end
end
