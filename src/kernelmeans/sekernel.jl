function calc_z(samples, prior, kernel)
    D = length(prior)
    z = samples .- Ref(mean(prior))
    Λ = _Λ(kernel)
    return scale(kernel) / sqrt(det(inv(Λ) * cov(prior) + I)) * exp.(-0.5 * PDMats.quad.(Ref(PDMat(inv(Λ + cov(prior)))), samples))
end

function calc_C(prior, kernel)
    D = length(prior)
    Λ = _Λ(kernel)
    return scale(kernel) / sqrt(det(2 * inv(Λ) * cov(prior) + I))
end


_Λ(kernel::TransformedKernel{<:SEKernel,<:ARDTransform}) = Diagonal(abs2.(inv.(kernel.transform.v)))
_Λ(kernel::TransformedKernel{<:SEKernel,<:ScaleTransform}) = inv(abs2(first(kernel.transform.s))) * I
_Λ(::SEKernel) = I