
prior(m::AbstractBQModel) = m.prior
logintegrand(m::AbstractBQModel) = m.logintegrand
integrand(m::AbstractBQModel) = x -> exp(logintegrand(m)(x))
logprior(m::AbstractBQModel) = x -> logpdf(prior(m), x)
logjoint(m::AbstractBQModel) = x -> logprior(m)(x) + logintegrand(m)(x)

"""
    BayesModel(prior, logintegrand)

Model inheriting from AbstractMCMC.AbstractModel.
`prior` should be a multivariate distribution from `Distributions.jl`
at the moment `prior` has to be a `MvNormal` but this will improved in a later version
`logintegrand` should be the log of the function to integrate.
"""
struct BayesModel{Tp,Ti} <: AbstractBQModel{Tp,Ti}
    prior::Tp
    logintegrand::Ti
end
