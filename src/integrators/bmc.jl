"""
    BMC(k::Kernel)

Bayesian Monte Carlo, the simplest approach for Bayesian Quadrature.
Assumes that the prior is Gaussian and that the integrand is as well.
Rasmussen, C. E. & Ghahramani, Z. Bayesian Monte Carlo. VDI Berichte 143–159 (2008).
"""
struct BMC{TK} <: AbstractIntegrator
    kernel::TK
end

function quadrature(bquad::BMC, model::AbstractBayesQuadModel{<:MvNormal}, samples)
    isempty(samples) && error("The collection of samples is empty")
    y = integrand(model).(samples)
    K = PDMat(kernelmatrix(kernel(bquad), samples))
    z = calc_z(samples, prior(model), kernel(bquad))
    C = calc_C(prior(model), kernel(bquad))
    return Normal(evaluate_mean(z, K, y), evaluate_var(z, K, C))
end

function evaluate_mean(z, K, y)
    return dot(z, K \ y)
end

function evaluate_var(z, K, C)
    return C - PDMats.invquad(K, z)
end