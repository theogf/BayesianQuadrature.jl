param(t::ScaleTransform) = inv(first(t.s))
param(t::ARDTransform) = inv(t.v)
param(t::LinearTransform) = inv(t.A)

include("sekernel.jl")