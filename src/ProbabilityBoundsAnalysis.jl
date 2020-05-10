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

using IntervalArithmetic, Distributions, Interpolations, Statistics, PyPlot, PyCall

import IntervalArithmetic: Interval, interval, AbstractInterval, isequal

import Base: show, -,
    +, *, /, //,
    <, >, ⊆, ^,
    rand

import Statistics: mean, var, std



#import Distributions: Normal, Beta, Uniform



## Global Variables (we may want to avoid these)

global steps = 200              # discretization levels of probability
global bOt = 0.001              # smallest quamtile to use if left tail is unbounded
global tOp = 0.999              # largest quamtile to use if right tail is unbounded
#global plotting = false         # if TRUE, p-pboxes are plotted whenever they are "show"
#global plottingEvery = false    # if TRUE, every p-box that's created is automatically plotted
#global cumulative = true        # if TRUE, plot CDFs, if FALSE, plot CCDFs (exceedance risks)
global pboxes = 0               # number of p-boxes that have been created
global verbose = 2              # how much warning messaging is wanted
global meanDisagreementMessage = "Disagreement between theoretical and observed mean\n"
global varDisagreementMessage = "Disagreement between theoretical and observed variance\n"

## Setter methods for globals

setSteps(t :: Int64) = (global steps = t)
setBOt(t :: Float64) = (global bOt = t)
setTOp(t:: Float64) = (global tOp = t)
#setPlotting(t :: Bool) = (global plotting = t)
#setPlottingEvery(t :: Bool) = (global plottingEvery = t)
#setCumulative(t :: Bool) = (global cumulative = t)
setPboxes(t :: Int64) = (global pboxes = t)
setVerbose(t:: Int) = (global verbose = t)
setMeanDisagreementMessage(t :: String) = (global meanDisagreementMessage = t)
setVarDisagreementMessage(t :: String) = (global varDisagreementMessage = t)


function showGlobals()

    println("\nShowing global variables:\n")
    println("Steps: $steps      ||  bOt: $bOt             ||  tOp: $tOp")
#    println("Plotting: $plotting  ||  plottingEvery: $plottingEvery   ||  cumulative: $cumulative")
    println("pboxes: $pboxes       ||  verbose: $verbose\n")
    print("meanDisagreementMessage: $meanDisagreementMessage")
    print("varDisagreementMessage:  $varDisagreementMessage")

end


export
    # Constructors
    Interval, interval, AbstractInterval, AbstractPbox, pbox, makepbox, 

    # Parametric
    normal, gaussian, N, Normal,
    U,uniform,
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

    minMaxMeanVar, MinMaxMeanVar, minmaxmeanvar,                            # Ferson
    Ferson, ferson, 

    # Moments
    mean, var, std, sideVariance, dwmean, dwVariance,

    # Arithmetic, binary
    conv, convFrechet, convPerfect, convOpposite, convFrechetNaive, balchProd,
    VKmeanproduct, VKmeanlo,VKmeanup,VKmeanlower, VKmeanupper,
    +, -, /, *,

    # Arithmetic, uniary
    negate, reciprocate, complement, shift,

    # Set based operations  
    env, imp, # Union, Intersection

    # Plots
    plot,

    # Checks
    ispbox, isinterval, isscalar, isvacuous, isequal,
    straddles, straddlingzero,

    # Sampling and cdf
    rand, cut, cdf,

    left, right, lefts, rights,

    # Copulas
    AbstractCopula, copula


include("pbox/pbox.jl")
include("pbox/dependency.jl")
include("pbox/arithmetic.jl")
include("pbox/MomentTransformations.jl")
include("pbox/distributions.jl")
include("pbox/special.jl")
include("pbox/plot_recipes.jl")
include("intervalStatistics/IntervalStatistics.jl")



end # module ProbabilityBoundsAnalysis