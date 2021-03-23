# ProbabilityBoundsAnalysis.jl
[![Build Status](https://travis-ci.com/AnderGray/ProbabilityBoundsAnalysis.jl.svg?branch=master)](https://travis-ci.com/github/AnderGray/ProbabilityBoundsAnalysis.jl)
[![codecov](https://codecov.io/gh/AnderGray/ProbabilityBoundsAnalysis.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/AnderGray/ProbabilityBoundsAnalysis.jl)
[![DOI](https://zenodo.org/badge/204744026.svg)](https://zenodo.org/badge/latestdoi/204744026)

[Probability bounds analysis](https://en.wikipedia.org/wiki/Probability_bounds_analysis) in Julia, a package for performing arithmetic between uncertain numbers. `ProbabilityBoundsAnalysis.jl` computes guaranteed bounds on functions of random variables, given only partial information about their marginals and dependence. Considered to be a form of rigorous computing with random variables.

<!---
This software began as a port from [pba.r](https://github.com/ScottFerson/pba.r) by Scott Ferson and Jason O'Rawe, Applied Biomathematics (2006)
--->

Supported uncertain numbers: 

  * scalars
  * [probability distributions](https://en.wikipedia.org/wiki/Cumulative_distribution_function)
  * [intervals](https://en.wikipedia.org/wiki/Interval_arithmetic)
  * [probability boxes](https://en.wikipedia.org/wiki/Probability_box) (p-boxes)

For more information and use, please see the [docs](https://andergray.github.io/ProbabilityBoundsAnalysis.jl/dev/)

Installation
---

`ProbabilityBoundsAnalysis.jl` is a registered Julia package, and so the latest release can be installed using the Julia package manager:

```julia
julia> ]
(v1.0) pkg> add ProbabilityBoundsAnalysis
julia> using ProbabilityBoundsAnalysis
```

### related packages:
* [pba.py](https://github.com/Institute-for-Risk-and-Uncertainty/pba-for-python): Python version of this software.
* [pba.r](https://github.com/ScottFerson/pba.r): R version of this software.
* [RAMAS® RiskCalc](https://www.ramas.com/riskcalc): a commerical software for distribution-free risk analysis.
* [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl): the interval arithmetic package used in this software.
* [ValidatedNumerics.jl](https://github.com/JuliaIntervals/ValidatedNumerics.jl): a suite of julia packages for rigorous computations.

Acknowledgements
---

The authors would like to thank the gracious support from the EPSRC iCase studentship award 15220067. We also acknowledge funding from UKRI via the EPSRC and ESRC Centre for Doctoral Training in Risk and Uncertainty Quantification and Management in Complex Systems. This research is funded by the Engineering Physical Sciences Research Council (EPSRC) with grant no. EP/R006768/1, which is greatly acknowledged for its funding and support. This work has been carried out within the framework of the EUROfusion Consortium and has received funding from the Euratom research and training programme 2014-2018 and 2019-2020 under grant agreement No 633053. The views and opinions expressed herein do not necessarily reflect those of the European Commission.

