@testset "kernelembeddings" begin 
    rng = MersenneTwister(42)
    N = 3
    measure = MvNormal(ones(3), ones(3))
    k = SqExponentialKernel()
    l = 2.0
    σ = 0.5
    ke = KernelEmbedding(k, σ, l, measure)
    kernel = σ * with_lengthscale(k, l)
    @test BQ.scale(ke) == σ

    sample = [rand(rng, 3)]
    @test kernel_mean(ke, sample) ≈ [mean(kernel.(sample, eachcol(rand(rng, measure, 10000))))] atol=1e-2
    @test kernel_variance(ke) ≈ mean(kernel.(eachcol(rand(rng, measure, 10000)), eachcol(rand(rng, measure, 10000)))) atol=1e-2

    s = 2.0
    l = [1.0, 2.0]
    L = LowerTriangular(rand(2, 2))
    @test BQ.Λ(s) ≈ s^2 * I
    @test BQ.Λ(l) ≈ Diagonal(abs2.(l))
    @test Matrix(BQ.Λ(L)) ≈ L * L'
end
