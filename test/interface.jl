@testset "Interface" begin
    rng = Random.MersenneTwister(42)
    D = 5
    prior = MvNormal(0.5 * ones(D))
    likelihood = MvNormal(7.0 * ones(D))
    integrand = x->pdf(likelihood, x)
    m = BayesModel(prior, integrand)
    s = PriorSampling()
    bquad = BayesQuad(transform(SEKernel(), 0.1))
    pZ, x  = bquad(m, s; nsamples=100)
    trueZ = exp(-0.5 * (logdet(cov(prior) + cov(likelihood)) + D * log(2π)))
    @test mean(pZ) ≈ trueZ atol=1e-3
end