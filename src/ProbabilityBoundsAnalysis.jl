######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Definition of the ProbabilityBoundsAnalysis.jl module
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
#
#   Port of R code pba.r by Scott Ferson and Jason O'Rawe, Applied Biomathematics
#   Origional code available at: https://github.com/ScottFerson/pba.r
######

#####
#       To install:
#
#			->  using package manager
#
#				julia> ]                                                     (opens package manager)
#               (v1.0) pkg> add ProbabilityBoundsAnalysis
#
#           ->  download package directly from github:
#
#               julia> ]                                                     (opens package manager)
#               (v1.0) pkg> add https://github.com/AnderGray/ProbabilityBoundsAnalysis.jl.git
#               julia> using ProbabilityBoundsAnalysis
#
#           ->  or if you have source code
#
#               julia> include("directory/of/source/src/ProbabilityBoundsAnalysis.jl")
#               julia> using Main.ProbabilityBoundsAnalysis
#
#
####

__precompile__(true)

module ProbabilityBoundsAnalysis

using Base: Float64
using Reexport, Distributions, Interpolations, PyCall, Distributed, PyPlot

@reexport using IntervalArithmetic, Distributed, Statistics, LinearAlgebra

#using PyPlot

import IntervalArithmetic: Interval, interval, AbstractInterval, isequal

import Base: show, -,
    +, *, /, //,
    <, >, ⊆, ^, intersect, issubset,
    rand, sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, exp, log, Threads.@spawn,
    Threads.@threads

import IntervalArithmetic: intersect, issubset

import Statistics: mean, var, std

import PyPlot: plot, scatter

import Distributions: cdf #Normal, Beta, Uniform


abstract type AbstractPbox <: Real end


## Global Variables (we may want to avoid these)

global steps = 200              # discretization levels of probability
global bOt = 0.001              # smallest quamtile to use if left tail is unbounded
global tOp = 0.999              # largest quamtile to use if right tail is unbounded
global verbose = 2              # how much warning messaging is wanted
global meanDisagreementMessage = "Disagreement between theoretical and observed mean\n"
global varDisagreementMessage = "Disagreement between theoretical and observed variance\n"

## Setter methods for globals

setSteps(t :: Int64) = (global steps = t)
setBOt(t :: Float64) = (global bOt = t)
setTOp(t:: Float64) = (global tOp = t)
setVerbose(t:: Int) = (global verbose = t)
setMeanDisagreementMessage(t :: String) = (global meanDisagreementMessage = t)
setVarDisagreementMessage(t :: String) = (global varDisagreementMessage = t)


function showGlobals()

    println("\nShowing global variables:\n")
    println("Steps: $steps      ||  bOt: $bOt             ||  tOp: $tOp")
    println("verbose: $verbose\n")
    print("meanDisagreementMessage: $meanDisagreementMessage")
    print("varDisagreementMessage:  $varDisagreementMessage")

end


export
    # Constructors
    Interval, interval, AbstractInterval, AbstractPbox, pbox, makepbox,

    # Parametric
    normal, gaussian, N,
    U, uniform,
    beta, betaPrime, biweight, cauchy,
    chi, chisq, cosine,epanechnikov, erlang,
    exponential, fDist, frechet, gamma, ksdist,
    laplace, levy, lognormal,

    # C-boxes
    KN,

    # Distribution free
    meanVar, MeanVar, minvar, meanVar, meanStd,                             # Chebyshev
    Chebyshev, chebyshev, cheb,

    meanMin, MeanMin, meanmin, meanMax, MeanMax, meanmax,                   # Markov
    Markov, markov,

    meanMinMax, MeanMinMax, meanminmax,                                     # Cantelli
    Cantelli, cantelli,

    minMeanVar, maxMeanVar,

    mmms, mmmv, minMaxMeanStd,

    minMaxMeanVar, MinMaxMeanVar, minmaxmeanvar,                            # Ferson
    Ferson, ferson,

    # Moments
    mean, var, std, sideVariance, dwmean, dwVariance,

    # Arithmetic, binary
    conv, convIndep, convFrechet, convPerfect, convOpposite, convFrechetNaive, balchProd,
    VKmeanproduct, VKmeanlo,VKmeanup,VKmeanlower, VKmeanupper, tauRho, sigma, convCorr,
    +, -, /, *,

    # Arithmetic, uniary
    negate, reciprocate, complement, shift,

    # Univariate functions
    sin, cos, tan, sinh, cosh, tanh, asin, acos, atan,
    exp, log, ^,

    # Set based operations
    env, imp,  intersect,# Union, Intersection

    # Plots
    plot,

    # Checks
    ispbox, isinterval, isscalar, isvacuous, isequal,
    straddles, straddlingzero,

    # Sampling and cdf
    rand, cut, cdf,

    left, right, lefts, rights,
    checkMomentsAndRanges,

    mixture,

    # Copulas
    AbstractCopula, copula, GauCopula, πCop, M, W, Frank, Clayton, Frechet, copulaSigma, copulaSigmaSlow, τCopula, KendalCopula, ρCopula, SpearmanCopula,
    copulaTau, copulaUnary, cdfU, cdfD, mass, conditionalX, conditionalY, bivpbox, plotStep, sample, sampleCond, scatter,
    plotBoxes


include("pbox/pbox.jl")
#include("pbox/dependency.jl")
include("pbox/copulas.jl")
include("pbox/NormalDistribution.jl")
include("pbox/arithmetic.jl")
include("pbox/MomentTransformations.jl")
include("pbox/distributions.jl")
include("pbox/special.jl")
include("pbox/plot_recipes.jl")
include("intervalStatistics/IntervalStatistics.jl")



end # module ProbabilityBoundsAnalysis
