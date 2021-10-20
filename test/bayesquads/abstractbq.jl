@testset "abstractbq" begin
    @test_throws ArgumentError BQ.check_kernel_parameters(-1.0, 2)
    @test_throws ArgumentError BQ.check_kernel_parameters([1.0, -1.0], 2)
    @test_throws ArgumentError BQ.check_kernel_parameters(rand(4, 4), 2)
    @test_throws ArgumentError BQ.check_kernel_parameters(1.0, -1)
    @test_nowarn BQ.check_kernel_parameters(1.0, 1.0)

    @test BQ.get_kernel_params(SqExponentialKernel(); a=1, b=2) == (SqExponentialKernel(), (a=1, b=2))
    @test BQ.get_kernel_params(3.0 * SqExponentialKernel(); l=1, σ=2) == (SqExponentialKernel(), (l=1, σ=3.0))
    @test BQ.get_kernel_params(3.0 * with_lengthscale(SqExponentialKernel(), 2.0); l=1, σ=2) == (SqExponentialKernel(), (l=2.0, σ=3.0))

    @test_throws ArgumentError BQ.check_transform(SelectTransform([1,3]))
    @test_nowarn BQ.check_transform(ScaleTransform(2.0))

    k = 3.0 * with_lengthscale(SqExponentialKernel(), 2.0)
    bquad = BayesQuad(k)
    x = rand(10)
    y = rand(10)
    @test BQ.kernel(bquad)(x, y) == k(x, y) # Kernels can not be compared directly
end