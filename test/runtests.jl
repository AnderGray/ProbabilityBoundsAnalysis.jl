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
include("pbox/comparisons.jl")
include("pbox/parameter_change.jl")
