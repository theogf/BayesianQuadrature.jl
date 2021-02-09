module BayesianQuadrature

using AbstractMCMC
using Distributions
using KernelFunctions
using LinearAlgebra
using PDMats
using Random

export bayesquad

abstract type AbstractIntegrator end
abstract type AbstractBayesQuadModel{Tp, Ti} <: AbstractMCMC.AbstractModel end

include(joinpath("integrators", "integrators.jl"))
include(joinpath("samplers", "samplers.jl"))
include(joinpath("kernelmeans", "kernels.jl"))
include("interface.jl")
include("models.jl")
include("utils.jl")

end

