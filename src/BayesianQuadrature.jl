module BayesianQuadrature

using AbstractMCMC
using Distributions
using KernelFunctions
using LinearAlgebra

export bayesquad

abstract type AbstractIntegrator end

include(joinpath("integrators", "integrators.jl"))
include(joinpath("samplers", "samplers.jl"))
include(joinpath("kernelmeans", "kernels.jl"))
include("utils.jl")

end

