# User guide

To run Bayesian quadrature multiple ingredients are needed:
- [The model](@ref), containing the integrand function `f` and the prior `p_0`
- [The sampling process](@ref), the algorithm choosing the samples ``\{x_i, y_i\}`` used to estimate the integral
- [The Bayesian quadrature algorithm](@ref), containing both the kernel and the way the posterior of the integral is computed 

## The model

So far only one model is implemented: the [`BayesModel`](@ref).
It takes a `MvNormal` prior and a log-likelihood `function`, in the future more options for the prior will be available using importance re-weighting!
In Bayesian terms, given the data ``y`` and latent parameters ``\theta``, with given likelihood ``p(y|\theta)`` and prior ``p_0(\theta)``,
we are interested in the evidence ``p(y) = \int p(y|\theta)p_0(\theta)d\theta``.
One needs to pass `p_0::MvNormal` ``\equiv p_0`` and `loglike(theta) = logpdf(likelihood(theta), y)` ``\equiv \log p(y|\theta)``   

## The sampling process

One need to select the different samples ``\{x_i,y_i\}``.
Since these samples do not need to belong to a specific distribution, all methods are allowed.
The options so far are:
- [`PriorSampling`](@ref), the most standard approach, sampling directly from `p_0`

## The Bayesian quadrature algorithm

This is the algorithm doing the heavy work. It describes the underlying Gaussian Process via 
a kernel from [KernelFunctions.jl](https://github.com/JuliaGaussianProcesses/KernelFunctions.jl).
When using a `SqExponentialKernel` with a `MvNormal` prior, the resulting distribution of the integral can be computed analytically.
The options are:
- [`BayesQuad`](@ref), this is the vanilla algorithm, which tries to compute the integral exactly.