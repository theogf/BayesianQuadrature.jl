function calc_z(samples, prior, bquad::BayesQuad{<:SqExponentialKernel})
    D = length(prior)
    z = samples .- Ref(mean(prior))
    B = Λ(bquad)
    return scale(bquad) / sqrt(det(inv(B) * cov(prior) + I)) *
           exp.(-0.5 * PDMats.quad.(Ref(PDMat(inv(B + cov(prior)))), samples))
end

function calc_C(prior, kernel)
    D = length(prior)
    B = Λ(bquad)
    return scale(bquad) / sqrt(det(2 * inv(B) * cov(prior) + I))
end