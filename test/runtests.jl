######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#  Running unit tests
#
######

using Test, ProbabilityBoundsAnalysis, IntervalArithmetic, Distributions, Documenter

include("pbox/pbox.jl")
include("pbox/distributions.jl")
include("pbox/arithmetic.jl")

DocMeta.setdocmeta!(ProbabilityBoundsAnalysis, :DocTestSetup, :(using ProbabilityBoundsAnalysis, IntervalArithmetic, Random; Random.seed!(31415)); recursive=true)
doctest(ProbabilityBoundsAnalysis, manual = false)
