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
using Distributions, Interpolations, Distributed, Requires

using IntervalArithmetic, Distributed, Statistics, LinearAlgebra

import IntervalArithmetic: Interval, interval, AbstractInterval, isequal

import Base: show, -,
    +, *, /, //,
    <, >, ⊆, ^, intersect, issubset,
    rand, sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, exp, log, Threads.@spawn,
    Threads.@threads

import IntervalArithmetic: intersect, issubset

import Statistics: mean, var, std

import Distributions: cdf #Normal, Beta, Uniform


abstract type AbstractPbox <: Real end


## Global Variables (we may want to avoid these)

mutable struct PboxParameters

    steps :: Int64              # discretization levels of probability
    bOt :: Float64              # smallest quamtile to use if left tail is unbounded
    tOp :: Float64              # largest quamtile to use if right tail is unbounded
    verbose :: Int64            # how much warning messaging is wanted

    PboxParameters() = new(200, 0.001, 0.999, 2)
end
## Setter methods for globals

const parametersPBA = PboxParameters()

setSteps(t :: Int64) = parametersPBA.steps = t
setBOt(t :: Float64) = parametersPBA.bOt = t
setTOp(t:: Float64) = parametersPBA.tOp = t
setVerbose(t:: Int) = parametersPBA.verbose = t


function show_parameters()

    println("\nShowing parameters:\n")
    println("Steps: $(parametersPBA.steps)      ||  bOt: $(parametersPBA.bOt)             ||  tOp: $(parametersPBA.tOp)")
    println("verbose: $(parametersPBA.verbose)\n")

end


export
    # Constructors
    Interval, interval, AbstractInterval, AbstractPbox, pbox, makepbox, parametersPBA,

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
include("pbox/copulas.jl")
include("pbox/NormalDistribution.jl")
include("pbox/arithmetic.jl")
include("pbox/MomentTransformations.jl")
include("pbox/distributions.jl")
include("pbox/special.jl")

function __init__()

    @require PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee" begin
        @require PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin

            using PyPlot, PyCall
            import PyPlot: plot, scatter

            pyimport_conda("mpl_toolkits.mplot3d", "mpl_toolkits")
            art3d = PyObject(PyPlot.art3D)
            mpl = pyimport("matplotlib")
            using3D()
            include("pbox/plots.jl")
        end
    end
end

end # module ProbabilityBoundsAnalysis
