"""
    PriorSampling()

Sampler which will use the prior distribution from the given model to provide samples.
"""
struct PriorSampling <: AbstractBQSampler end

function AbstractMCMC.step(rng::Random.AbstractRNG, model, ::PriorSampling; kwargs...)
    return rand(rng, p_0(model)), nothing
end

function AbstractMCMC.step(
    rng::Random.AbstractRNG, model, ::PriorSampling, ::Any; kwargs...
)
    return rand(rng, p_0(model)), nothing
end
