using BayesianQuadrature
using KernelFunctions
using Distributions
using Statistics
D = 2
p_0 = MvNormal(ones(D)) # As for now the prior must be a MvNormal
likelihood = MvNormal(0.5 * ones(D))
log_f(x) = logpdf(likelihood, x) # The logarithm of the Integrand log_f, the log-likelihood function typically
model = BayesModel(p_0, log_f) # Combine both to create the model
logbquad = LogBayesQuad(SEKernel(); l=0.1, σ=1.0) # Will simply approximate p(y|x) with a GP (only works with SEKernel for now
bquad = BayesQuad(SEKernel(); l=0.1, σ=1.0) # Will simply approximate p(y|x) with a GP (only works with SEKernel for now
sampler = PriorSampling() # Will sample from the prior p_0
n_s = 100
log_p_I, _ = logbquad(model, sampler; nsamples=n_s) # Returns a Normal distribution
p_I, _ = bquad(model, sampler; nsamples=n_s+100) # Returns a Normal distribution
@show p_I # Normal{Float64}(μ=0.07063602778449946, σ=0.0028050929209120458)
trueZ = exp(-0.5 * (logdet(cov(p_0) + cov(likelihood)) + D * log(2π)))
mean(p_I)
mean(log_p_I)