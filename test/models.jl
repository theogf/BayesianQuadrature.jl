@testset "Models logfunctions" begin
    rng = Random.MersenneTwister(42)
    D = 5
    prior = MvNormal(0.5 * ones(D))
    likelihood = MvNormal(7.0 * ones(D))
    integrand = x->pdf(likelihood, x)
    m = BayesModel(prior, integrand)  

    x1 = randn(rng,Float64,D)
    x2 = randn(rng,Float64,D)
    x3 = randn(rng,Float64,D)
    x4 = randn(rng,Float64,D)
    @test logprior(m)(x1) ≈ logpdf(prior,x1)
    @test logprior(m)(x2) ≈ logpdf(prior,x2)
    @test logprior(m)(x3) ≈ logpdf(prior,x3)
    @test logprior(m)(x4) ≈ logpdf(prior,x4)
    @test logjoint(m)(x1) ≈ logpdf(prior,x1) + logpdf(likelihood,x1)
    @test logjoint(m)(x2) ≈ logpdf(prior,x2) + logpdf(likelihood,x2)
    @test logjoint(m)(x3) ≈ logpdf(prior,x3) + logpdf(likelihood,x3)
    @test logjoint(m)(x4) ≈ logpdf(prior,x4) + logpdf(likelihood,x4)
end