using BayesianQuadrature
using Distributions
using KernelFunctions
using LinearAlgebra
using Random
using Test
const BQ = BayesianQuadrature

@testset "BayesianQuadrature.jl" begin
    include("interface.jl")

    @info "Testing Integrators"
    @testset "Integrators" begin
        include(joinpath("integrators", "integrators.jl"))
        include(joinpath("integrators", "bmc.jl"))
    end

    @info "Testing Samplers"
    @testset "Samplers" begin
        include(joinpath("samplers", "samplers.jl"))
        include(joinpath("samplers", "montecarlo.jl"))
    end

    @info "Testing kernel means"
    @testset "Kernel Means" begin
        include(joinpath("kernelmeans", "kernels.jl"))
        include(joinpath("kernelmeans", "sekernel.jl"))
    end
end
