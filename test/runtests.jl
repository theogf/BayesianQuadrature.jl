using BayesianQuadrature
using Distributions
using KernelFunctions
using LinearAlgebra
using Random
using Test
const BQ = BayesianQuadrature

include("testing_tools.jl")

@testset "BayesianQuadrature.jl" begin
    include("interface.jl")
    include("models.jl")

    @info "Testing BayesQuads"
    @testset "BayesQuads" begin
        include(joinpath("bayesquads", "abstractbq.jl"))
        include(joinpath("bayesquads", "bayesquad.jl"))
        include(joinpath("bayesquads", "logbayesquad.jl"))
    end

    @info "Testing Samplers"
    @testset "Samplers" begin
        include(joinpath("samplers", "abstractbqsampler.jl"))
        include(joinpath("samplers", "priorsampling.jl"))
    end

    @info "Testing kernel means"
    @testset "Kernel Means" begin
        include(joinpath("kernelmeans", "kernels.jl"))
        include(joinpath("kernelmeans", "sekernel.jl"))
    end
end
