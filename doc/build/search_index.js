var documenterSearchIndex = {"docs":
[{"location":"#ProbabilityBoundsAnalysis.jl-Documentation","page":"Home","title":"ProbabilityBoundsAnalysis.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"mutable struct pbox <: AbstractPbox\n\npbox(u::Union{Missing, Array{<:Real}, Real} = missing, d=u; shape = \"\", name = \"\", ml= missing, mh = missing, vl = missing, vh = missing, interpolation = \"linear\",bob = missing, perfect = missing, opposite = missing, depends = missing, dids=missing, bounded = [false, false])\n\npbox( x :: Array{Interval{T}, 1}, bounded = [true, true]) where T <: Real\npbox( x :: Array{T, 2}, bounded = [true, true]) where T <: Real\n\ncdf(s :: pbox, x :: Real)\ncdf(s :: pbox, x::Interval)\n\nmass(s :: pbox, lo :: Real, hi :: Real)\nmass(s :: pbox, x:: Interval)\n\nmakepbox(x...)\n\nmean(x::pbox)\n\nvar(x::pbox)\n\nstd(x::pbox)","category":"page"}]
}