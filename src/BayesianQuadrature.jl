module BayesianQuadrature

using AbstractGPs
using Distributions
using KernelFunctions

export bayesquad

abstract type AbstractSelectiveSampler end

abstract type AbstractIntegrator end

end

