@testset "bayesquad" begin
    k = SqExponentialKernel()

    σ = 4.0
    @test BQ.scale(BayesQuad(σ * k)) == σ

    test_convergence(BayesQuad(SqExponentialKernel()), PriorSampling())
end
