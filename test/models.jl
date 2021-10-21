@testset "Models logfunctions" begin
    rng = Random.MersenneTwister(42)
    D = 5
    p0 = MvNormal(0.5 * ones(D))
    likelihood = MvNormal(7.0 * ones(D))
    logf = x -> logpdf(likelihood, x)
    m = BayesModel(p0, logf)

    x1 = randn(rng, Float64, D)
    @test integrand(m)(x1) ≈ pdf(likelihood, x1)
    @test logprior(m)(x1) ≈ logpdf(p0, x1)
    @test logjoint(m)(x1) ≈ logpdf(p0, x1) + logpdf(likelihood, x1)
    @test BQ.p_0(m) == p0
    @test logintegrand(m) == logf
end
