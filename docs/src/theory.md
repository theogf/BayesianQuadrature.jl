# Theory

!!! note
    This is a very quick and dirty introduction to Bayesian quadrature. For a more concrete introduction,
    see these [nice slides](https://www.cs.toronto.edu/~duvenaud/talks/intro_bq.pdf)

Bayesian Quadrature is a probabilistic numerics methods to compute integrals.
The idea is to compute the integral
```math
I = \int f(x) p(x)dx,
```
where ``p(x)\sim \mathcal{N}``, by representing the function ``f(x)`` by a [Gaussian Process](https://en.wikipedia.org/wiki/Gaussian_process).

However, to know what the posterior of ``f`` is like, we need samples ``\{x_i, y_i\}`` where ``y_i = f(x_i) + \epsilon``, with ``epsilon`` being some noise.
In the end, instead of looking for an approximation of ``I``, we obtain a posterior distribution for it values:

```math
    p(I|\{x_i,y_i\}) = \mathcal{N}(\mu, \sigma^2)
```

Under some conditions, the mean ``\mu`` and the variance ``\sigma^2`` can be found analytically.