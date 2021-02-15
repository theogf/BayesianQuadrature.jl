function (bquad::AbstractBayesQuad)(
    rng::Random.AbstractRNG,
    model::AbstractBayesQuadModel,
    sampler::AbstractMCMC.AbstractSampler;
    x_init = [],
    nsamples = 200,
)
    x = sample(rng, model, sampler, nsamples)
    return quadrature(bquad, model, x), x
end


function (bquad::AbstractBayesQuad)(
    model::AbstractBayesQuadModel,
    sampler::AbstractMCMC.AbstractSampler;
    x_init = [],
    nsamples = 200,
)
    bquad(GLOBAL_RNG, model, sampler; x_init=x_init, nsamples=nsamples)
end