#=
    Compute z, where z_i = ∫k(x, x_i)p(x)dx, i.e. the kernel mean also
    called the kernel embedding of a distribution
=#
function calc_z(samples, prior, bquad::AbstractBQ{<:SqExponentialKernel})
    z = samples .- Ref(mean(prior))
    B = Λ(bquad.l)
    return scale(bquad) / sqrt(det(inv(B) * cov(prior) + I)) *
           exp.(-0.5 * PDMats.invquad.(Ref(PDMat(B + cov(prior))), samples))
end

#=
    Compute C = ∫∫k(x, x')p(x)p(x')dxdx', i.e. the kernel variance
=#
function calc_C(prior, bquad::AbstractBQ{<:SqExponentialKernel})
    B = Λ(bquad.l)
    return scale(bquad) / sqrt(det(2 * inv(B) * cov(prior) + I))
end
