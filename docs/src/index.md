# ProbabilityBoundsAnalysis.jl


[Probability bounds analysis](https://en.wikipedia.org/wiki/Probability_bounds_analysis) in Julia, a package for performing arithmetic between uncertain numbers. `ProbabilityBoundsAnalysis.jl` computes guaranteed bounds on functions of random variables, given only partial information about their marginals and dependence. Considered to be a form of rigorous computing with random variables.


Supported uncertain numbers: 

  * scalars
  * [probability distributions](https://en.wikipedia.org/wiki/Cumulative_distribution_function)
  * [intervals](https://en.wikipedia.org/wiki/Interval_arithmetic)
  * [probability boxes](https://en.wikipedia.org/wiki/Probability_box) (p-boxes)

---

### Authors

- [Ander Gray](https://www.researchgate.net/profile/Ander_Gray), Institute for Risk and Uncertainty, University of Liverpool
- [Scott Ferson](https://www.liverpool.ac.uk/engineering/staff/scott-ferson/), Institute for Risk and Uncertainty, University of Liverpool


### Collaborators

- [Marco De Angelis](https://marcodeangelis.github.io), Institute for Risk and Uncertainty, University of Liverpool
- [Nick Gray](https://riskinstitute.uk/people/nickgray/), Institute for Risk and Uncertainty, University of Liverpool
- [Alexander Wimbush](https://riskinstitute.uk/people/alexanderwimbush/), Institute for Risk and Uncertainty, University of Liverpool

---

Installation
---
Two ways to install and use:

**1. From the julia package manager**

You may download the lastest release by:
```julia
julia> ]
(v1.0) pkg> add ProbabilityBoundsAnalysis
julia> using ProbabilityBoundsAnalysis
```

or the latest version by:

```julia
julia> ]
(v1.0) pkg> add ProbabilityBoundsAnalysis#master
julia> using ProbabilityBoundsAnalysis
```

or

```julia
julia> ]
(v1.0) pkg> add https://github.com/AnderGray/ProbabilityBoundsAnalysis.jl.git
julia> using ProbabilityBoundsAnalysis
```

**2. Downloading the source code**
```julia
git clone https://github.com/AnderGray/ProbabilityBoundsAnalysis.jl.git

julia> include("ProbabilityBoundsAnalysis.jl/src/ProbabilityBoundsAnalysis.jl")
julia> using Main.ProbabilityBoundsAnalysis
```

---

### related packages:
* [pba.py](https://github.com/Institute-for-Risk-and-Uncertainty/pba-for-python): Python version of this software.
* [pba.r](https://github.com/ScottFerson/pba.r): R version of this software.
* [RAMAS® RiskCalc](https://www.ramas.com/riskcalc): a commerical software for distribution-free risk analysis.
* [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl): the interval arithmetic package used in this software.
* [ValidatedNumerics.jl](https://github.com/JuliaIntervals/ValidatedNumerics.jl): a suite of julia packages for rigorous computations.

---

### Acknowledgements

The authors would like to thank the gracious support from the EPSRC iCase studentship award 15220067. We also gratefully acknowledge funding from UKRI via the EPSRC and ESRC Centre for Doctoral Training in Risk and Uncertainty Quantification and Management in Complex Systems. This research is funded by the Engineering Physical Sciences Research Council (EPSRC) with grant no. EP/R006768/1, which is greatly acknowledged for its funding and support. This work has been carried out within the framework of the EUROfusion Consortium and has received funding from the Euratom research and training programme 2014-2018 and 2019-2020 under grant agreement No 633053. The views and opinions expressed herein do not necessarily reflect those of the European Commission.
