
prior(m::AbstractBayesQuadModel) = m.prior
integrand(m::AbstractBayesQuadModel) = m.integrand
logprior(m::AbstractBayesQuadModel) = x->logpdf(prior(m), x)
logjoint(m::AbstractBayesQuadModel) = x->logprior(m)(x) + log(integrand(m)(x))

"""
    BayesModel(prior, integrand)

Model inheriting from AbstractMCMC.AbstractModel.
`prior` should be a multivariate distribution from `Distributions.jl`
at the moment `prior` has to be a `MvNormal` but this will improved in a later version
`integrand` should be the function to integrate.
"""
struct BayesModel{Tp,Ti} <: AbstractBayesQuadModel{Tp,Ti}
    prior::Tp
    integrand::Ti
end
