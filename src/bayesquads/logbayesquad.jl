"""
    LogBayesQuad(k::Kernel; l=1.0, σ::Real=1.0)

Bayesian Quadrature object.
You can pass any kernel and the lengthscale and variance will be extracted.
`l` can be a `Real`, a `AbstractVector` or a `LowerTriangular`.
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
    x_c = sample_candidates(bquad, samples)
    logf = logintegrand(model).(samples)
    f = exp.(logf)

    f_c_0 = mean.(predict(bquad, samples, f, x_c))
    logf_c_0 = mean.(predict(bquad, samples, logf, x_c))
    Δ_c = exp.(logf_c_0) - f_c_0
    
    z = calc_z(samples, prior(model), bquad)
    K = kernelpdmat(kernel(bquad), samples)
    C = calc_C(prior(model), bquad)
    
    z_c = calc_z(x_c, prior(model), bquad)
    K_c = kernelpdmat(kernel(bquad), x_c)
    m_evidence = evaluate_mean(z, K, f)
    m_correction = evaluate_mean(z_c, K_c, Δ_c)

    var_evidence = evaluate_var(z, K, C)
    var_correction = evaluate_var(z_c, K_c, C)
    return Normal(m_evidence + m_correction, sqrt(var_evidence + var_correction))
end

function predict(bquad::LogBayesQuad, samples, f, x_c)
    gp = GP(kernel(bquad))
    fx = gp(samples)
    f_post = posterior(fx, f)
    return f_c = mean(f_post(x_c)) 
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