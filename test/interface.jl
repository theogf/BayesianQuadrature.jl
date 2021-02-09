# @testset "Interface" begin
rng = Random.MersenneTwister(42)
D = 5
prior = MvNormal(ones(D))
likelihood = MvNormal(0.1 * ones(D))
integrand = x->pdf(likelihood, x)
m = BQ.BayesModel(prior, integrand)
s = BQ.MonteCarlo()
i = BQ.BMC(transform(SEKernel(), 10.0))
pZ, x  = BQ.bayesquad(m, i, s; nsamples=100)
Σ = inv(inv(cov(prior)) + inv(cov(likelihood)))
posterior = MvNormal(Diagonal(Σ))
trueZ = 1 / sqrt(det(2π * Σ))
@test mean(bqZ) ≈ trueZ atol=1e-1