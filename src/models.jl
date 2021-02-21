
prior(m::AbstractBayesQuadModel) = m.prior
logintegrand(m::AbstractBayesQuadModel) = m.logintegrand
integrand(m::AbstractBayesQuadModel) = x->exp(logintegrand(m)(x))
logprior(m::AbstractBayesQuadModel) = x->logpdf(prior(m), x)
logjoint(m::AbstractBayesQuadModel) = x->logprior(m)(x) + logintegrand(m)(x)

"""
    BayesModel(prior, logintegrand)

Model inheriting from AbstractMCMC.AbstractModel.
`prior` should be a multivariate distribution from `Distributions.jl`
at the moment `prior` has to be a `MvNormal` but this will improved in a later version
`logintegrand` should be the logarithm function to integrate.
"""
struct BayesModel{Tp,Ti} <: AbstractBayesQuadModel{Tp,Ti}
    prior::Tp
    logintegrand::Ti
end
