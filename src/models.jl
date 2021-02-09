struct BayesModel{Tp,Ti} <: AbstractBayesQuadModel{Tp,Ti}
    prior::Tp
    integrand::Ti
end

prior(m::AbstractBayesQuadModel) = m.prior
integrand(m::AbstractBayesQuadModel) = m.integrand
logprior(m::AbstractBayesQuadModel) = x->logpdf(prior(m), x)
logjoint(m::AbstractBayesQuadModel) = x->logprior(m)(x) + log(integrand(x))