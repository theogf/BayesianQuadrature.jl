struct MonteCarlo <: AbstractMCMC.AbstractSampler end

function AbstractMCMC.step(rng::Random.AbstractRNG, model, ::MonteCarlo; kwargs...)
    return rand(rng, prior(model)), nothing
end

function AbstractMCMC.step(rng::Random.AbstractRNG, model, ::MonteCarlo, ::Any; kwargs...)
    return rand(rng, prior(model)), nothing
end