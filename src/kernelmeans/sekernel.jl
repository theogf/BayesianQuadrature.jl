function calc_z(samples, prior, kernel)
    D = length(prior)
    z = samples .- Ref(mean(prior))
    Λ = _Λ(kernel)
    return scale(kernel) / sqrt(det(inv(Λ) * cov(prior) + I)) * exp.(-0.5 * quad.(Ref(inv(Λ + cov(prior)))))
end

function calc_C(prior, kernel)
    D = length(prior)
    Λ = _Λ(kernel)
    return scale(kernel) / sqrt(det(2 * inv(Λ) * cov(prior) + I))
end


_Λ(kernel::TransformedKernel{<:SEKernel,<:ARDTransform}) = Diagonal(inv.(kernel.transform.v))
_Λ(kernel::TransformedKernel{<:SEKernel,<:ScaleTransform}) = inv(first(kernel.transform.s)) * I
_Λ(::SEKernel) = I