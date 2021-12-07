module BayesianQuadrature

using AbstractGPs
using AbstractMCMC
using Distributions
using KernelFunctions
using LinearAlgebra
using PDMats
using Random

export BayesQuad, LogBayesQuad
export PriorSampling
export BayesModel
export prior, integrand, logintegrand, logprior, logjoint
export BQ # Short version for calling BayesianQuadrature
export KernelEmbedding
export kernel_mean, kernel_variance

const BQ = BayesianQuadrature

"""
    AbstractBQ

General class of models to perform quadrature.
It needs to implement `quadrature(::AbstractBQ, model::AbstractBQModel, samples)`
"""
abstract type AbstractBQ{TK} end

"""
    AbstractBQSampler

General class for sampling.
Should implement the `AbstracMCMC` API
"""
abstract type AbstractBQSampler <: AbstractMCMC.AbstractSampler end

"""
    AbstractBQModel

General model class.
Should implement `p_0(m)` and `logintegrand(m)`
"""
abstract type AbstractBQModel{Tp,Ti} <: AbstractMCMC.AbstractModel end

include("bayesquads/abstractbq.jl")
include("samplers/abstractbqsampler.jl")
include("kernelembeddings/kernelembedding.jl")
include("interface.jl")
include("models.jl")
include("utils.jl")

end
