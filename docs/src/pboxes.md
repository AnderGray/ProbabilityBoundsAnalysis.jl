# Constructing probability boxes


A probability distribution can be created by using it's shape and parameters:

```julia
julia> a = normal(0,1)
Pbox: 	  ~ normal ( range=[-3.090232,3.090232], mean=0.0, var=1.0)
```
In `ProbabilityBoundsAnalysis.jl` probability distributions are p-boxes, but with equal bounds and precise moments.

The `IntervalArithmetic.jl` package is used to create intervals:

```julia
julia> b = interval(0,1)
[0, 1]
```

`IntervalArithmetic.jl` will work natively with `ProbabilityBoundsAnalysis.jl`, however an interval can be converted to a p-box in the following way:

```julia
julia> b = makepbox(interval(0,1))
Pbox: 	  ~  ( range=[0.0,1.0], mean=[0.0,1.0], var=[0.0,0.25])
```

There are a number of ways that p-boxes can be created. For example they are the result of arithemtic between random variables with unknown dependence. They can also be defined by using a distributions shape but with interval parameters:

```julia
julia> c = normal(interval(0,1),1)
Pbox: 	  ~ normal ( range=[-3.09023,4.0902322], mean=[0.0,1.0], var=1.0)
```
or by taking the envelope over a number of uncertain numbers:
```julia
julia> d = normal(-1,1); 
julia> e = normal(1, interval(1,2));
julia> f = env(d,e)
Pbox: 	  ~ normal ( range=[-5.18046,7.1804646], mean=[-1.0,1.0], var=[1.0,4.0])
```
and may be plotted as follows:

```julia
julia> plot(f)
```
![alt text](https://github.com/AnderGray/ProbabilityBoundsAnalysis.jl/tree/master/docs/plots/PbaPlot1.png "a probability box")

In `ProbabilityBoundsAnalysis.jl` all plots of uncertain numbers are of their cdfs.

Supported parametric distributions:

<!-----

|           |         |        |              |         |             |
|:---------:|:-------:|:------:|:------------:|:-------:|:-----------:|
|   normal  | uniform |  beta  |   betaPrime  | biweght |    Cauchy   |
|    chi    |  chisq  | cosine | epanechnikov |  erlang | exponential |
|   fDist   | frechet |  gamma |    ksdist    | laplace |     levy    |
| lognormal |         |        |              |         |             |

Supported distribution free p-boxes:

|           |         |        |              |         |
|:---------:|:-------:|:------:|:------------:|:-------:|
|   meanVar  | meanMin |  meanMax  |   meanMinMax  | minMaxMeanVar|

----->
* normal
* uniform
* beta
* betaPrime
* biweght 
* cauchy
* chi
* chisq
* cosine
* epanechnikov 
* erlang
* exponential

Supported distribution free p-boxes:

* meanVar
* meanMin
* meanMax
* meanMinMax
* minMaxMeanVar

KN c-boxes also supported.

All constructors support interval arguments.

