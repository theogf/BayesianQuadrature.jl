module BayesianQuadrature

using AbstractGPs
using Distributions
using KernelFunctions

export bayesquad

abstract type AbstractSelectiveSampler end

abstract type AbstractIntegrator end

include(joinpath("integrators", "integrators.jl"))
include(joinpath("samplers", "samplers.jl"))
include(joinpath("kernelmeans", "kernels.jl"))
include("utils.jl")

end

