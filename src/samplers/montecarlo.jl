struct MonteCarlo{Tp} <: AbstractSampler
    p::Tp
end

function step(sampler::MonteCarlo, ::AbstractVector, ::AbstractVector, args...)
    return rand(sampler.p)
end