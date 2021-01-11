# Modelling dependencies with copulas

Any possible multivariate dependence can be encoded in a [copula](https://en.wikipedia.org/wiki/Copula_(probability_theory)). Copulas, sometimes called dependency functions, are a joint cdf with standard uniform marginals, and are a way to model dependency independently of marginal distributions. 

The supported copulas are:

  * π: independence
  * M: perfect dependence and the upper Frechét bound
  * W: opposite dependence and the lower Frechét bound
  * Gaussian: family with correlation coefficient r (-1 for W, 0 for π, 1 for M)
  * Frank: family with parameter s (0 for W, 1 for π, Inf for M)
  * Clayton: family with parameter t (-1 for W, 0 for π, Inf for M)

A copula can be created and plotted in the following way:

```julia
julia> C = Gaussian(0.7)
Copula ~ Gau(r=0.7)

julia> plotCdf(a)
```
![alt text](https://github.com/AnderGray/ProbabilityBoundsAnalysis.jl/blob/master/doc/plots/GaussianCopula.png "Gaussian copula with a correlation of 0.7")

A copula C is the function C:[0,1]<sup>d</sup> &rarr;[0,1], where d is the dimension of the copula. Only bivariate copulas are considered so far in `ProbabilityBoundsAnalysis.jl`. A copula can be evaluated and sampled in the following way:


```julia
julia> C(0.5,0.5)
0-dimensional Array{Float64,0}:
0.37340834444668247

julia> C([0.2,0.6],[0.3,0.7])
2×2 Array{Float64,2}:
 0.143315  0.194937
 0.273398  0.52667 
 
julia> samps = sample(C,10^6);
```

Given any marginal distributions F<sub>X</sub>(x) and F<sub>Y</sub>(y) and a copla C( : ), a joint distribution can be created: H(x,y) = C(F<sub>X</sub>(x),F<sub>Y</sub>(y)).  In `ProbabilityBoundsAnalysis.jl` a joint distribution can be created by passing marginals to the copula. For example, a distribution with beta marginals and a gaussian copula:

```julia
julia> J = C(Beta(4,2),Beta(4,2))
Joint ~ Gau( r=0.7, Beta{Float64}(α=4.0, β=2.0), Beta{Float64}(α=4.0, β=2.0) )

julia> plotDensity(J)
```
![alt text](https://github.com/AnderGray/ProbabilityBoundsAnalysis.jl/blob/master/doc/plots/JointDist1.png "Bivariate Beta with gaussian copula with correlation of 0.7")

A joint distribution can also be sampled:

```julia
julia> Jsamps = sample(J,10^6)
julia> scatter(Jsamps)
```
![alt text](https://github.com/AnderGray/ProbabilityBoundsAnalysis.jl/blob/master/doc/plots/Jsamples.png "Scatter plot of samples of J")