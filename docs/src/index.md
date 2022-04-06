# ProbabilityBoundsAnalysis.jl


[Probability bounds analysis](https://en.wikipedia.org/wiki/Probability_bounds_analysis) in Julia, a package for performing arithmetic between uncertain numbers. `ProbabilityBoundsAnalysis.jl` computes guaranteed bounds on functions of random variables, given only partial information about their marginals and dependence. Considered to be a form of rigorous computing with random variables.


Supported uncertain numbers: 

  * scalars
  * [probability distributions](https://en.wikipedia.org/wiki/Cumulative_distribution_function)
  * [intervals](https://en.wikipedia.org/wiki/Interval_arithmetic)
  * [probability boxes](https://en.wikipedia.org/wiki/Probability_box) (p-boxes)

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

For plotting
---

[PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl) is used for all of the plotting functionality, it is however not installed by default. If you wish to plot p-boxes, you will need to install and the specifically call `PyPlot`.

**To install**
```julia
julia> ]
(v1.0) pkg> add PyPlot
```

`Matplotlib` needs to be installed on your computer. If it is not, it can be installed by

```julia
julia -e 'ENV["PYTHON"]=""; using Pkg; Pkg.add("Conda"); using Conda; Conda.add("matplotlib"); Pkg.add("PyCall"); Pkg.build("PyCall"); Pkg.add("PyPlot");'
```
**To use**

```julia
julia> using ProbabilityBoundsAnalysis, PyPlot
julia> a = normal(0,1)
julia> plot(a)
```

---

Community Guidelines
---
### Contributing
If you like this package or find it useful, do consider contributing. Any improvements or suggestions for new features are welcome, either to the package or this documentation. If you wish to make a contribution, please     make a push changes to a branch with a descriptive name e.g,  `agray-new_feature` and make a pull request. Please write a description of your modifications/additions in the pull request.

### Unit testing

If you make a change to the code base, consider also writing a test. The tests can be found in `test/runtests.jl`, which are run in every push to a pull request using github's CI. If you wish to run the tests locally
```julia
julia --project --color=yes test/runtests.jl
```
### Reporting issues and seeking support
Issues and bugs can be reported using the `issues` tab on the repository. Feel free to report any issues that you have. Questions about functionality, and new feature requests may also be communicated through `issues`.

Further support can be forwarded to `ander.gray@ukaea.uk`.

---

### related packages:
* [pba.py](https://github.com/Institute-for-Risk-and-Uncertainty/pba-for-python): Python version of this software.
* [pba.r](https://github.com/ScottFerson/pba.r): R version of this software.
* [RAMAS® RiskCalc](https://www.ramas.com/riskcalc): a commerical software for distribution-free risk analysis.
* [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl): the interval arithmetic package used in this software.
* [ValidatedNumerics.jl](https://github.com/JuliaIntervals/ValidatedNumerics.jl): a suite of julia packages for rigorous computations.

---

### Acknowledgements

Financial support from the EPSRC iCase studentship award 15220067 is gratefully acknowledged. We also gratefully acknowledge funding from UKRI via the EPSRC and ESRC Centre for Doctoral Training in Risk and Uncertainty Quantification and Management in Complex Systems. This research is funded by the Engineering Physical Sciences Research Council (EPSRC) with grant no. EP/R006768/1, which is greatly acknowledged for its funding and support. This work has been carried out within the framework of the EUROfusion Consortium and has received funding from the Euratom research and training programme 2014-2018 and 2019-2020 under grant agreement No 633053. The views and opinions expressed herein do not necessarily reflect those of the European Commission.
