"""
    PriorSampling()

Sampler which will use the prior distribution from the given model to provide samples.
"""
struct PriorSampling <: AbstractMCMC.AbstractSampler end

function AbstractMCMC.step(rng::Random.AbstractRNG, model, ::PriorSampling; kwargs...)
    return rand(rng, prior(model)), nothing
end

function AbstractMCMC.step(
    rng::Random.AbstractRNG, model, ::PriorSampling, ::Any; kwargs...
)
    return rand(rng, prior(model)), nothing
end
