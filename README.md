# pba.jl
[Probability bounds analysis](https://en.wikipedia.org/wiki/Probability_bounds_analysis) in Julia, a package for performing arithmetic between uncertain numbers. `pba.jl` computes guaranteed bounds on functions of random variables, given only partial information about their marginals and dependence. Considered to be a form of rigorous computing with random variables.

The mayority of this code is a port from [pba.r](https://github.com/AnderGray/pba.jl) by Scott Ferson and Jason O'Rawe, Applied Biomathematics

Supported uncertain numbers: 

  * scalars
  * probability distributions
  * [intervals](https://en.wikipedia.org/wiki/Interval_arithmetic)
  * [probability boxes](https://en.wikipedia.org/wiki/Probability_box)

Supported arithmetic operations:


Supported dependent arithmetic between uncertain numbers:

|                           |     independent    | dependency known   | dependency unknown | perfect/opposite     | partial information  |
|---------------------------|:------------------:|--------------------|--------------------|----------------------|----------------------|
| intervals                 | not known to exist | no solution exists |         yes        |  solutions exist |  solutions exist |
| probability distributions |         yes        |         tbc        |         yes        |          yes         |          solutions exist         |
| probability boxes         |         yes        |         tbc        |         yes        |          yes         |          solutions exist         |

**tbc: to be continued. General solutions exist and will be implemented.**


Installation
---
Two ways to install and use:

1. From the julia package manager
```julia
julia> ]
(v1.0) pkg> add https://github.com/AnderGray/pba.jl.git
julia> using pba
```
(will be `add pba` when full package is released)

2. Downloading Source code
```julia
julia> include("directory/of/source/src/pba.jl")
julia> using Main.pba
```

Once installed, uncertain numbers can be created

Uncertain numbers
---

A probability distribution can be using it's shape and parameters:

```julia
julia> a = normal(0,1)
Pbox: 	  ~ normal ( range=[-3.090232306167818,3.090232306167818], mean=0.0, var=1.0)
```
In `pba.jl` probability distributions are p-boxes, but with equal bounds and precise moments.

The `IntervalArithmetic.jl` package is used to create intervals:

```julia
julia> b = interval(0,1)
[0, 1]
```

`IntervalArithmetic.jl` will work natively with `pba.jl`, however an interval can be converted to a p-box in the following way:

```julia
julia> b = makepbox(interval(0,1))
Pbox: 	  ~  ( range=[0.0,1.0], mean=[0.0,1.0], var=[0.0,0.25])
```

There are a number of ways that p-boxes can be created. For example they are the result of arithemtic between random variables with unknown dependence. They can also be defined by using a distributions shape but with interval parameters:

```julia
julia> c = normal(interval(0,1),1)
Pbox: 	  ~ normal ( range=[-3.090232306167818,4.0902323061678185], mean=[0.0,1.0], var=1.0)
```
or by taking the envelope over a number of uncertain numbers:
```julia
julia> d = normal(-1,1); 
julia> e = normal(1, interval(1,2));
julia> f = env(d,e)
Pbox: 	  ~ normal ( range=[-5.180464612335636,7.180464612335636], mean=[-1.0,1.0], var=[1.0,4.0])
```
and may be plotted as follows:

```julia
julia> plotpbox(f)
```


### Arithmetic


Interval statistics 
---

For performing descriptive statistics of interval data sets. Not yet a part of the main module

To use:
```julia
julia> include("directory/of/source/src/IntervalStatistics/IntervalStatistics.jl")
```


Dependency modelling
----



Bibliography
---

For probability bound analysis:
* [*Ferson, S., V. Kreinovich, L. Ginzburg, K. Sentz and D.S. Myers. 2003. Constructing probability boxes and Dempster-Shafer structures. Sandia National Laboratories, SAND2002-4015, Albuquerque, New Mexico*](https://www.osti.gov/servlets/purl/1427258)
* *Ferson, S., and J. Siegrist. 2012. Verified computation with probabilities. Pages 95-122 in Uncertainty Quantification in Scientific Computing, edited by A. Dienstfrey and R.F. Boisvert, Springer, New York*
* [*Beer, M., S. Ferson, and V. Kreinovich. 2013. Imprecise probabilities in engineering analyses. Mechanical Systems and Signal Processing 37: 429*](https://digitalcommons.utep.edu/cgi/viewcontent.cgi?article=1733&=&context=cs_techrep&=&sei-redir=1&referer=https%253A%252F%252Fscholar.google.com%252Fscholar%253Fhl%253Den%2526as_sdt%253D0%25252C5%2526q%253DBeer%25252C%252BM.%25252C%252BS.%252BFerson%25252C%252Band%252BV.%252BKreinovich.%252B2013.%252BImprecise%252Bprobabilities%252Bin%252Bengineering%252Banalyses.%252BMechanical%252BSystems%252Band%252BSignal%252BProcessing%252B37%25253A%252B429%2526btnG%253D#search=%22Beer%2C%20M.%2C%20S.%20Ferson%2C%20V.%20Kreinovich.%202013.%20Imprecise%20probabilities%20engineering%20analyses.%20Mechanical%20Systems%20Signal%20Processing%2037%3A%20429%22)
* *Ferson, S., D.R.J. Moore, P. van der Brink, T.L. Estes, K. Gallagher, R. O’Connor and F. Verdonck. 2010. Bounding uncertainty analyses. Application of Uncertainty Analysis to Ecological Risks of Pesticides*

For interval statistics:
* [*Ferson, Scott, et al. "Experimental uncertainty estimation and statistics for data having interval uncertainty." Sandia National Laboratories, Report SAND2007-0939 162 (2007)*](https://www.researchgate.net/file.PostFileLoader.html?id=52b1b418d3df3e110f8b45b1&assetKey=AS%3A272184528834564%401441905256417)

### related packages:
* [puffin](https://github.com/AnderGray/pba.jl): (needs link) uncertainty compiler. Converts scientific codes in to one readable by a package like this.
* [pba.r](https://github.com/AnderGray/pba.jl): (needs link) R version of this software.
* [RAMAS® RiskCalc](https://www.ramas.com/riskcalc): a comerical software for distribution-free risk analysis.
* [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl): the interval arithmetic package used in this software.
* [ValidatedNumerics.jl](https://github.com/JuliaIntervals/ValidatedNumerics.jl): a suite of julia packages for rigorous computations.

Acknowledgements
---

For all those who believe in honest computation (add funders)
