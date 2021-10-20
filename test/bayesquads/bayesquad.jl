@testset "bayesquad" begin
    s = 2.0
    l = [1.0, 2.0]
    L = LowerTriangular(rand(2, 2))
    k = SqExponentialKernel()
    @test BQ.Λ(BayesQuad(k ∘ ScaleTransform(s))) ≈ 1/s^2 * I
    @test BQ.Λ(BayesQuad(k ∘ ARDTransform(l))) ≈ Diagonal(inv.(abs2.(l)))
    @test BQ.Λ(BayesQuad(k ∘ LinearTransform(L))) ≈ inv(L) * inv(L)'

    σ = 4.0
    @test BQ.scale(BayesQuad(σ * k)) == σ

    test_convergence(BayesQuad(SqExponentialKernel()), PriorSampling())
end