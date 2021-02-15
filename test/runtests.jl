using BayesianQuadrature
using Distributions
using KernelFunctions
using LinearAlgebra
using Random
using Test
const BQ = BayesianQuadrature

@testset "BayesianQuadrature.jl" begin
    include("interface.jl")

    @info "Testing BayesQuads"
    @testset "BayesQuads" begin
        include(joinpath("bayesquads", "bayesquads.jl"))
        include(joinpath("bayesquads", "bayesquad.jl"))
    end

    @info "Testing Samplers"
    @testset "Samplers" begin
        include(joinpath("samplers", "samplers.jl"))
        include(joinpath("samplers", "priorsampling.jl"))
    end

    @info "Testing kernel means"
    @testset "Kernel Means" begin
        include(joinpath("kernelmeans", "kernels.jl"))
        include(joinpath("kernelmeans", "sekernel.jl"))
    end
end
