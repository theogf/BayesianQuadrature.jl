function test_convergence(bquad, sampler; D=3)
    rng = Random.MersenneTwister(42)
    prior = MvNormal(0.5 * ones(D))
    likelihood = MvNormal(7.0 * ones(D))
    logintegrand = x->logpdf(likelihood, x)
    m = BayesModel(prior, logintegrand)
    pZ, x  = bquad(m, sampler; nsamples=100)
    trueZ = exp(-0.5 * (logdet(cov(prior) + cov(likelihood)) + D * log(2π)))
    @test mean(pZ) ≈ trueZ rtol=1e-2
end