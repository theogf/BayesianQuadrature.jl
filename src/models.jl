#= 
    General helping functions
=#
logprior(m::AbstractBQModel) = Base.Fix1(logpdf, p_0(m))
prior(m::AbstractBQModel) = exp ∘ logprior(m)

integrand(m::AbstractBQModel) = exp ∘ logintegrand(m)
logjoint(m::AbstractBQModel) = x -> logprior(m)(x) + logintegrand(m)(x)

"""
    BayesModel(prior, logintegrand) <: AbstractBQModel

Model inheriting from `AbstractMCMC.AbstractModel`.
`prior` should be a multivariate distribution from `Distributions.jl`
at the moment `prior` has to be a `MvNormal` but this will improved in a later version
`logintegrand` should be the log of the function to integrate.
"""
struct BayesModel{Tp,Ti} <: AbstractBQModel{Tp,Ti}
    prior::Tp
    logintegrand::Ti
end

p_0(m::BayesModel{<:AbstractMvNormal}) = m.prior
p_0(m::BayesModel) = MvNormal(ones(dim(m.prior)))
logintegrand(m::BayesModel{<:AbstractMvNormal}) = m.logintegrand
function logintegrand(m::AbstractBQModel)
    return function reweightedlogintegrand(x)
        return m.logintegrand(x) + logpdf(m.prior, x) - logprior(m)(x)
    end
end
