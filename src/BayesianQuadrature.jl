module BayesianQuadrature

using AbstractMCMC
using Distributions
using KernelFunctions
using LinearAlgebra
using PDMats
using Random

export BayesQuad, PriorSampling, BayesModel
export prior, integrand, logintegrand, logprior, logjoint
export BQ # Short version for calling BayesianQuadrature

const BQ = BayesianQuadrature

abstract type AbstractBQ end
abstract type AbstractBQSampler <: AbstractMCMC.AbstractSampler end
abstract type AbstractBQModel{Tp,Ti} <: AbstractMCMC.AbstractModel end

include(joinpath("bayesquads", "bayesquads.jl"))
include(joinpath("samplers", "samplers.jl"))
include(joinpath("kernelmeans", "kernels.jl"))
include("interface.jl")
include("models.jl")
include("utils.jl")

end
