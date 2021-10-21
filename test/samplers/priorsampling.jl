@testset "priorsampling" begin 
    sampler = PriorSampling()
    model = BayesModel(MvNormal(ones(2)), x->sum(x))
    samples = sample(model, sampler, 10000)
    @test mean(samples) ≈ zeros(2) atol=1e-1
    @test cov(samples) ≈ Matrix(I(2)) atol=1e-1
end
