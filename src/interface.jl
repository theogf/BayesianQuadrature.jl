"""


"""
function bayesquad(
    rng::Random.AbstractRNG,
    model::AbstractBayesQuadModel,
    integrator::AbstractIntegrator,
    sampler::AbstractSampler,
    x_init::AbstractVector;
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
    return quadrature(integrator, integrand, prior, samples)
end