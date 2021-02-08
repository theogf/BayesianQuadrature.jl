_Λ(::T) where {T<:Kernel} = error("Method _Λ not defined for kernels of type $T")
scale(kernel::ScaledKernel) = first(kernel.σ²)
scale(kernel::Kernel) = 1
scale(kernel::TransformedKernel{<:ScaledKernel}) = scale(kernel.kernel)

include("sekernel.jl")