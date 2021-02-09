"""


"""
function bayesquad(
    rng::Random.AbstractRNG,
    model::AbstractBayesQuadModel,
    integrator::AbstractIntegrator,
    sampler::AbstractMCMC.AbstractSampler;
    x_init = [],
    nsamples = 200,
)
    # logπ = logjoint(model)
    x = sample(rng, model, sampler, nsamples)
    # T = typeof(step(sampler, x_init, [], logπ))
    # samples = Vector{T}(undef, n_samples + 1)
    # samples[1] = x_init
    # for i in 1:n_samples
    #     samples[i+1] = step(sampler, samples[i], samples, logπ)
    # end
    return quadrature(integrator, model, x), x
end


function bayesquad(
    model::AbstractBayesQuadModel,
    integrator::AbstractIntegrator,
    sampler::AbstractMCMC.AbstractSampler;
    x_init = [],
    nsamples = 200,
)
    bayesquad(GLOBAL_RNG, model, integrator, sampler; x_init=x_init, nsamples=nsamples)
end