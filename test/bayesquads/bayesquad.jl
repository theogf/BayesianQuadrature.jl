@testset "bayesquad" begin
    s = 2.0
    l = [1.0, 2.0]
    L = LowerTriangular(rand(2, 2))
    k = SqExponentialKernel()
    @test BQ.Λ(s) ≈ s^2 * I
    @test BQ.Λ(l) ≈ Diagonal(abs2.(l))
    @test BQ.Λ(L) ≈ L * L'

    σ = 4.0
    @test BQ.scale(BayesQuad(σ * k)) == σ

    test_convergence(BayesQuad(SqExponentialKernel()), PriorSampling())
end