function calc_z(samples, prior, bquad::AbstractBayesQuad{<:SqExponentialKernel})
    D = length(prior)
    z = samples .- Ref(mean(prior))
    B = Λ(bquad)
    return scale(bquad) / sqrt(det(inv(B) * cov(prior) + I)) *
           exp.(-0.5 * PDMats.invquad.(Ref(PDMat(B + cov(prior))), samples))
end

function calc_C(prior, bquad::AbstractBayesQuad{<:SqExponentialKernel})
    D = length(prior)
    B = Λ(bquad)
    return scale(bquad) / sqrt(det(2 * inv(B) * cov(prior) + I))
end
