# BayesianQuadrature

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://theogf.github.io/BayesianQuadrature.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://theogf.github.io/BayesianQuadrature.jl/dev)
[![Build Status](https://github.com/theogf/BayesianQuadrature.jl/workflows/CI/badge.svg)](https://github.com/theogf/BayesianQuadrature.jl/actions)
[![Coverage](https://codecov.io/gh/theogf/BayesianQuadrature.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/theogf/BayesianQuadrature.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

Package for implementing different methods of Bayesian quadrature.
Bayesian quadrature consists in estimating the integral `I = ∫ f(x) p(x) dx` by using Gaussian Processes where `p(x)` is assumed to be Gaussian.
More precisely we replace `f(x)` by a GP by estimating `f` for multiple samples `x_i`.
We then get a posterior distribution for the integral : `p(I|{x_i}) = N(m, S)`.

Given a Bayesian problem `p(x|y) = p(y|x) p_0(x) / p(y)` you can estimate `p(y)` by calling :

```julia
prior = MvNormal(ones(2)) # As for now the prior must be a MvNormal
integrand(x) = pdf(MvNormal(0.5 * ones(2)), x) # Integrand, the likelihood function typically
model = BayesModel(prior, integrand)
integrator = BMC(SEKernel()) # Will simply approximate p(y|x) with a GP (only works with SEKernel for now
sampler = MCSampling() # Will sample from the prior p_0
I = bayesquad(model, integrator, sampler; nsamples=100) # Returns a Normal distribution
```
