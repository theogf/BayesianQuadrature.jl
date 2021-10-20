
prior(m::AbstractBQModel) = m.prior
logintegrand(m::AbstractBQModel) = m.logintegrand
integrand(m::AbstractBQModel) = x -> exp(logintegrand(m)(x))
logprior(m::AbstractBQModel) = x -> logpdf(prior(m), x)
logjoint(m::AbstractBQModel) = x -> logprior(m)(x) + logintegrand(m)(x)

"""
    BayesModel(prior, logintegrand) <: AbstractBayesQuadModel

Model inheriting from `AbstractMCMC.AbstractModel`.
`prior` should be a multivariate distribution from `Distributions.jl`
at the moment `prior` has to be a `MvNormal` but this will improved in a later version
`logintegrand` should be the log of the function to integrate.
"""
struct BayesModel{Tp,Ti} <: AbstractBQModel{Tp,Ti}
    prior::Tp
    logintegrand::Ti
end

_prior(m::AbstractBayesQuadModel) = m.prior
priord(m::AbstractBayesQuadModel{<:MvNormal}) = _prior(m)
priord(m::AbstractBayesQuadModel) = MvNormal(ones(length(_prior(m))))
prior(m::AbstractBayesQuadModel) = exp ∘ logprior(m)
logprior(m::AbstractBayesQuadModel) = x->logpdf(priord(m), x)

_logintegrand(m::AbstractBayesQuadModel) = m.logintegrand
logintegrand(m::AbstractBayesQuadModel{<:MvNormal}) = _logintegrand(m)
function logintegrand(m::AbstractBayesQuadModel)
    return function reweightedlogintegrand(x)
        _logintegrand(m)(x) + logpdf(_prior(m), x) - logprior(m)(x)
    end
end
integrand(m::AbstractBayesQuadModel) = exp ∘ logintegrand(m)

logjoint(m::AbstractBayesQuadModel) = x->logprior(m)(x) + logintegrand(m)(x)
