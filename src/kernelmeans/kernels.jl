_Λ(::T) where {T<:Kernel} = error("Method _Λ not defined for kernels of type $T")
scale(kernel::ScaleTransform) = first(kernel.σ²)
scale(kernel::Kernel) = 1
scale(::T) where {T<:Kernel} = error("Method scale not defined for kernels of type $T")

include("sekernel.jl")