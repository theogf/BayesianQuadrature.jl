@testset "Interface" begin
    rng = Random.MersenneTwister(42)
    D = 5
    prior = MvNormal(0.5 * ones(D))
    likelihood = MvNormal(7.0 * ones(D))
    integrand = x->pdf(likelihood, x)
    m = BQ.BayesModel(prior, integrand)
    s = BQ.MonteCarlo()
    i = BQ.BMC(transform(SEKernel(), 0.1))
    pZ, x  = BQ.bayesquad(m, i, s; nsamples=100)
    Σ = cov(prior) * inv(cov(prior) + cov(likelihood)) * cov(likelihood)
    posterior = MvNormal(Diagonal(Σ))
    trueZ = exp(-0.5 * (logdet(cov(prior) + cov(likelihood)) + D * log(2π)))
    @test mean(pZ) ≈ trueZ atol=1e-3
end