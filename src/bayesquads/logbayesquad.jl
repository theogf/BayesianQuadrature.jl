"""
    LogBayesQuad(k::Kernel; n_candidates=100, l=1.0, σ::Real=1.0)
Tool for running Bayesian Quadrature by assuming a GP on the log integrand instead
of the integrand.
`n_candidates` will be sampled around the first samples to perform a linear approximation
of the result.
See : 
- **Active Learning of Model Evidence Using Bayesian Quadrature** - Osborne et al. - NIPS 2012.
"""
struct LogBayesQuad{TK,Tl,Tσ} <: AbstractBayesQuad{TK, Tl}
    kernel::TK
    n_candidates::Int
    l::Tl
    σ::Tσ
end

function LogBayesQuad(k::Kernel; n_candidates=100, l=1.0, σ::Real=1.0)
    σ > 0 || ArgumentError("σ should be positive")
    l isa AbstractMatrix && (
        l isa LowerTriangular || throw(
            ArgumentError(
                "For l an AbstractMatrix, only LowerTriangular matrices are accepted"
            ),
        )
    )
    return LogBayesQuad(k, n_candidates, l, σ)
end

function LogBayesQuad(k::TransformedKernel{TK,Tt}; l=nothing, σ=1.0) where {TK,Tt}
    Tt <: Union{ScaleTransform,ARDTransform,LinearTransform} ||
        error("No lengthscale could be extracted from kernel $k,\n
        only ScaleTransform, ARDTransform and LinearTransform are allowed")
    l = param(k.transform)

    return LogBayesQuad(k.kernel; l=l, σ=σ)
end

function LogBayesQuad(k::ScaledKernel; l=1.0, σ=nothing)
    σ = first(k.σ²)
    return LogBayesQuad(k.kernel; l=l, σ=σ)
end

function quadrature(
    bquad::LogBayesQuad{<:SqExponentialKernel},
    model::AbstractBayesQuadModel{<:MvNormal},
    samples,
)
    isempty(samples) && error("The collection of samples is empty")
    nsamples = length(samples)
    logf = logintegrand(model).(samples) # Evaluate integrand on samples
    normalize = false
    if normalize
        max_logf = maximum(logf) # Find maximum
        logf .-= max_logf # Normalize log integrand
    end
    f = exp.(logf) # Compute integrand

    x_c = sample_candidates(bquad, samples) # Sample candidates around the samples

    gp = create_gp(bquad, samples)
    f_c_0 = mean.(predict(gp, f, x_c)) # Predict integrand on x_c
    logf_c_0 = mean.(predict(gp, logf, x_c)) # Predict log-integrand on x_c
    Δ_c = exp.(logf_c_0) - f_c_0 # Compute difference of predictions
    
    C = calc_C(prior(model), bquad)

    z = calc_z(samples, prior(model), bquad)
    K = kernelpdmat(kernel(bquad), samples)
    
    z_c = calc_z(x_c, prior(model), bquad)
    K_c = kernelpdmat(kernel(bquad), x_c)

    m_evidence = evaluate_mean(z, K, f)
    m_correction = evaluate_mean(z_c, K_c, Δ_c)

    var_evidence = evaluate_var(z, K, C)
    var_correction = evaluate_var(z_c, K_c, C)
    if normalize
        m = exp(log(m_evidence + m_correction) + max_logf)
        v = exp(log(var_evidence + var_correction) + 2 * max_logf)
    else
        m = m_evidence + m_correction
        v = var_evidence + var_correction
    end
    return Normal(m, sqrt(v))
end

function predict(gp::AbstractGPs.FiniteGP, f, x_c)
    f_post = posterior(gp, f)
    return f_c = mean(f_post(x_c)) 
end

function create_gp(bquad::LogBayesQuad, samples)
    gp = GP(kernel(bquad))
    fx = gp(samples)    
end

function sample_candidates(bquad::LogBayesQuad, samples)
    n_c = bquad.n_candidates
    n_samples = length(samples)
    n_dim = length(first(samples))
    directions = [rand(n_dim) .- 0.5 |> x -> x / norm(x, 1) for _ in 1:n_c]
    directions = lengthscale(bquad) .* directions
    p_indices = rand(1:n_samples, n_c)
    return directions .+= samples[p_indices] 
end