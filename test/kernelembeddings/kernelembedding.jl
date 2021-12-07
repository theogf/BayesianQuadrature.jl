@testset "kernelembeddings" begin 
    rng = MersenneTwister(42)
    N = 3
    measure = MvNormal(ones(3), ones(3))
    k = SqExponentialKernel()
    l = 2.0
    σ = 0.5
    ke = KernelEmbedding(k, l, σ, measure)
    kernel = σ * with_lengthscale(k, l)
    @test scale(ke) == σ

    sample = [rand(rng, 3)]
    @test kernel_mean(ke, sample) ≈ [mean(kernel.(sample, eachcol(rand(rng, measure, 10000))))]
end
