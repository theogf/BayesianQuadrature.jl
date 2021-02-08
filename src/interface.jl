"""


"""
function bayesquad(
    integrand,
    prior,
    integrator::AbstractIntegrator,
    sampler::AbstractSampler,
    x_init::AbstractVector;
    n_samples = 200,
)
    logπ(x) = logpdf(prior, x) + log(integrand(x))
    T = typeof(step(sampler, x_init, [], logπ))
    samples = Vector{T}(undef, n_samples + 1)
    samples[1] = x_init
    for i in 1:n_samples
        samples[i+1] = step(sampler, samples[i], samples, logπ)
    end
    return quadrature(integrator, integrand, prior, samples)
end