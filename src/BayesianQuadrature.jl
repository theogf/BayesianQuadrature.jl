module BayesianQuadrature

using AbstractGPs
using AbstractMCMC
using Distributions
using KernelFunctions
using LinearAlgebra
using PDMats
using Random

export BayesQuad, LogBayesQuad
export PriorSampling, BayesModel
export prior, integrand, logintegrand, logprior, logjoint

abstract type AbstractBayesQuad{Kernel, T} end
abstract type AbstractBayesQuadModel{Tp, Ti} <: AbstractMCMC.AbstractModel end

include(joinpath("bayesquads", "bayesquads.jl"))
include(joinpath("samplers", "samplers.jl"))
include(joinpath("kernelmeans", "kernels.jl"))
include("interface.jl")
include("models.jl")
include("utils.jl")

end

