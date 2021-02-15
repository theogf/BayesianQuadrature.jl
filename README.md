# BayesianQuadrature

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://theogf.github.io/BayesianQuadrature.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://theogf.github.io/BayesianQuadrature.jl/dev)
[![Build Status](https://github.com/theogf/BayesianQuadrature.jl/workflows/CI/badge.svg)](https://github.com/theogf/BayesianQuadrature.jl/actions)
[![Coverage Status](https://coveralls.io/repos/github/theogf/BayesianQuadrature.jl/badge.svg?branch=main)](https://coveralls.io/github/theogf/BayesianQuadrature.jl?branch=main)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

Package for implementing different methods of Bayesian quadrature.
Bayesian quadrature consists in estimating the integral `I = ∫ f(x) p(x) dx` by using Gaussian Processes where `p(x)` is assumed to be Gaussian.
More precisely we replace `f(x)` by a GP by estimating `f` for multiple samples `x_i`.
We then get a posterior distribution for the integral : `p(I|{x_i}) = N(m, S)`.

Given a Bayesian problem `p(x|y) = p(y|x) p_0(x) / p(y)` you can estimate `p(y)` by calling :

```julia
using BayesianQuadrature
using Distributions
using KernelFunctions
prior = MvNormal(ones(2)) # As for now the prior must be a MvNormal
integrand(x) = pdf(MvNormal(0.5 * ones(2)), x) # Integrand, the likelihood function typically
model = BayesModel(prior, integrand) # Combine both to create the model
bquad = BayesQuad(SEKernel(); l=0.1, σ=1.0) # Will simply approximate p(y|x) with a GP (only works with SEKernel for now
sampler = PriorSampling() # Will sample from the prior p_0
p_I, _ = bquad(model, sampler; nsamples=100) # Returns a Normal distribution
@show p_I # Normal{Float64}(μ=0.07063602778449946, σ=0.0028050929209120458)
```
