# API

```@meta
CurrentModule = ProbabilityBoundsAnalysis
DocTestSetup = quote
    using ProbabilityBoundsAnalysis
end
```

```@docs

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

uniform(min :: Interval, max :: Interval)

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