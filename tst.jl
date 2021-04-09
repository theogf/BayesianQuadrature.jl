using BayesianQuadrature
using Distributions
using KernelFunctions
using Random

rng = Random.MersenneTwister(42)

p_0 = MvNormal(ones(2)) # As for now the prior must be a MvNormal
log_f(x) = logpdf(MvNormal(0.5 * ones(2)), x) # The logarithm of the Integrand log_f, the log-likelihood function typically
m = BayesModel(p_0, log_f) # Combine both to create the model
bquad = BayesQuad(SEKernel(); l=10.0, σ=1.0) # Will simply approximate p(y|x) with a GP (only works with SEKernel for now
sampler = PriorSampling() # Will sample from the prior p_0
p_I, _ = bquad(rng, m, sampler; nsamples=200) # Returns a Normal distribution
@show p_I 