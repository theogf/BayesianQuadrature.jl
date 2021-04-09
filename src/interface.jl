function (bquad::AbstractBQ)(
    rng::Random.AbstractRNG,
    model::AbstractBQModel,
    sampler::AbstractMCMC.AbstractSampler;
    x_init=[],
    nsamples=200,
    callback=nothing,
)
    x = sample(rng, model, sampler, nsamples; callback=callback, bquad=bquad)
    return quadrature(bquad, model, x), x
end

function (bquad::AbstractBQ)(
    model::AbstractBQModel,
    sampler::AbstractMCMC.AbstractSampler;
    x_init=[],
    nsamples=200,
    callback=nothing,
)
    return bquad(
        GLOBAL_RNG, model, sampler; x_init=x_init, nsamples=nsamples, callback=callback
    )
end
