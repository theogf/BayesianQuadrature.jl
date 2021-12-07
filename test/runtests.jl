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

    @info "Testing kernel embeddings"
    @testset "Kernel Embeddings" begin
        include(joinpath("kernelembeddings", "kernelembedding.jl"))
    end
end
