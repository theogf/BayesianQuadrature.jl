using BayesianQuadrature
using KernelFunctions
using Distributions
p_0 = MvNormal(ones(2)) # As for now the prior must be a MvNormal
log_f(x) = logpdf(MvNormal(0.5 * ones(2)), x) # The logarithm of the Integrand log_f, the log-likelihood function typically
model = BayesModel(p_0, log_f) # Combine both to create the model
bquad = LogBayesQuad(SEKernel(); l=0.1, σ=1.0) # Will simply approximate p(y|x) with a GP (only works with SEKernel for now
sampler = PriorSampling() # Will sample from the prior p_0
p_I, _ = bquad(model, sampler; nsamples=100) # Returns a Normal distribution
@show p_I # Normal{Float64}(μ=0.07063602778449946, σ=0.0028050929209120458)