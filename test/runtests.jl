######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#  Running unit tests
#
######

using Test, ProbabilityBoundsAnalysis, IntervalArithmetic, Distributions, Documenter

include("pbox/pbox.jl")
include("pbox/special.jl")
include("pbox/comparisons.jl")
include("pbox/distributions.jl")
include("pbox/copulas.jl")
include("pbox/arithmetic.jl")
include("pbox/parameter_change.jl")
include("pbox/plots.jl")
