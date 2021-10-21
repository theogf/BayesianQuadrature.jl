@testset "logbayesquad" begin
    test_convergence(LogBayesQuad(SqExponentialKernel()), PriorSampling())
end
