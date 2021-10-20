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
export BQ # Short version for calling BayesianQuadrature

const BQ = BayesianQuadrature

"""
    AbstractBQ

General class of models to perform quadrature.
It needs to implement `quadrature(::AbstractBQ, model::AbstractBQModel, samples)`
"""
abstract type AbstractBQ end
abstract type AbstractBQSampler <: AbstractMCMC.AbstractSampler end
abstract type AbstractBQModel{Tp,Ti} <: AbstractMCMC.AbstractModel end

include("bayesquads/abstractbq.jl")
include("samplers/abstractbqsampler.jl")
include("kernelmeans/kernels.jl")
include("interface.jl")
include("models.jl")
include("utils.jl")

end
