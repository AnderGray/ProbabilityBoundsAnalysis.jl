# API


```@docs
pbox(u::Union{Missing, Array{<:Real}, Real} = missing, d=u; shape = "", name = "", ml= missing, mh = missing, vl = missing, vh = missing, interpolation = "linear",bob = missing, perfect = missing, opposite = missing, depends = missing, dids=missing, bounded = [false, false])

pbox( x :: Array{Interval{T}, 1}, bounded = [true, true]) where T <: Real
pbox( x :: Array{T, 2}, bounded = [true, true]) where T <: Real

cdf(s :: pbox, x :: Real)
cdf(s :: pbox, x::Interval)

mass(s :: pbox, lo :: Real, hi :: Real)
mass(s :: pbox, x:: Interval)

makepbox(x...)

mean(x::pbox)

var(x::pbox)

std(x::pbox)

env(x...)

imp(x...)

normal(mean, std)

uniform(Min, Max)

beta(α, β)

lognormal(  μ, θ)

KN(k, n)

cut(x, p :: Real)


```

## Index

```@index
Pages = ["api.md"]
Module = ["ProbabilityBoundsAnalysis"]
```