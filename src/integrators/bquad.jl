"""



Rasmussen, C. E. & Ghahramani, Z. Bayesian Monte Carlo. VDI Berichte 143–159 (2008).
"""
struct BQuad{TGP} <: AbstractIntegrator
    gp::TGP
end

function BQuad(k::Kernel, μ₀::MeanFunction=ZeroMean())
    BQuad(GP(k, μ₀))
end

function quadrature(bquad::BQuad, integrand, prior, samples)
    y = integrand.(samples)
    K = cov(gp(bquad), samples)
    z = calc_z(x, prior, kernel(bquad))
    C = calc_C(prior, kernel(bquad))
    return Normal(int_mean(z, K, y), int_variance(z, K, C))
end