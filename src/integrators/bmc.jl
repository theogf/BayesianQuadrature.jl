"""
    BMC(k::Kernel)

Bayesian Monte Carlo, the simplest approach for Bayesian Quadrature.
Assumes that the prior is Gaussian and that the integrand is as well.
Rasmussen, C. E. & Ghahramani, Z. Bayesian Monte Carlo. VDI Berichte 143–159 (2008).
"""
struct BMC{TK} <: AbstractIntegrator
    kernel::TK
end

function quadrature(bquad::BMC, integrand, prior::MvNormal, samples)
    isempty(samples) && error("The collection of samples is empty")
    y = integrand.(samples)
    K = cov(gp(bquad), samples)
    z = calc_z(x, prior, kernel(bquad))
    C = calc_C(prior, kernel(bquad))
    return Normal(int_mean(z, K, y), int_variance(z, K, C))
end

function evaluate_mean(z, K, y)
    return dot(z, K \ y)
end

function evaluate_var(z, K, C)
    return C - invquad(K, z)
end