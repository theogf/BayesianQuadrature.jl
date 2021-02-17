@testset "Models logfunctions" begin
    rng = Random.MersenneTwister(42)
    D = 5
    p0 = MvNormal(0.5 * ones(D))
    likelihood = MvNormal(7.0 * ones(D))
    integrand = x->pdf(likelihood, x)
    m = BayesModel(prior, integrand)  

    x1 = randn(rng,Float64,D)
    @test logprior(m)(x1) ≈ logpdf(prior,x1)
    @test logjoint(m)(x1) ≈ logpdf(prior,x1) + logpdf(likelihood,x1)
    @test BQ.prior(m) == prior
    @test BQ.integrand(m) == integrand
end
