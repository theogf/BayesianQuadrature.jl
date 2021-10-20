@testset "bayesquad" begin
    @testset "lengthscale and variance" begin
        σ²_k = 1.5; l_k = 1.5
        k = σ²_k * SqExponentialKernel() ∘ ScaleTransform(1. /  l_k)
        bq = BayesQuad(k)
        @test bq.σ ≈ σ²_k
        @test bq.l ≈ l_k
    end
end