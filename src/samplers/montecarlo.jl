struct MonteCarlo{Tp} <: AbstractMCMC.AbstractSampler end

function AbstractMCMC.step(rng::Random.AbstractRNG, model, ::MonteCarlo; kwargs...)
    return rand(rng, prior(model))
end

function AbstractMCMC.step(rng::Random.AbstractRNG, model, ::MonteCarlo, ::Any; kwargs...)
    return rand(rng, prior(model))
end