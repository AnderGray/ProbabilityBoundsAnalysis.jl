######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Definition of bivariate copula class, joint distribution class and methods
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
######

###
#   Notes: This is a script for performing tests with copulas for plotting, sampling and
#   finding resulting copulas from arithemtic between uncertain variables.
#
#   Would it be possible to represent the set of all copulas as an expansion of orthogonal
#   polynomials? Could be used as a sentivitivty analysis for dependancy. Something equivilant
#   to disperive Monte Carlo, where not only the correlatin coefficient is changed but the entire
#   dependancy structure. Could this expansion be used in the arithmetic?
#
#   Or can we do it as sums of the Frechet bounds? (Shuffles of M)
#   See: Nelson page 68.
#
#   For algebra of random variables, including algebra of covaraince See https://en.wikipedia.org/wiki/Algebra_of_random_variables
#       -> The above link is also useful for moment propagation. Could you perform covariance algebra when the covarainces are partially known
#
#
#   One model for dependency trancking, is to assume a copula family, and compute the covariance of the results of opartaions with their inputs
#   using the arithemtic of covaraince described above. For this we would need an interval covaraince algebra, and a way of performing convolutions with it.
#
###

#=
using Distributions
#using PyPlot
using PyCall
using LinearAlgebra
using Interpolations
include("NormalDistribution.jl")
=#
#global  n  = 200;        # Number of descretizations
#global  pn = 200;        # Number of descretization for plotting
#global bOt = 0.001;              # smallest quamtile to use if left tail is unbounded
#global tOp = 0.999;              # largest quamtile to use if right tail is unbounded

abstract type AbstractCopula <: Real end
abstract type AbstractJoint <: Real end

mutable struct copula <: AbstractCopula

    cdfU ::  Array{Float64,2}
    cdfD ::  Array{Float64,2}
    family
    param

    function copula(cdfU = missing, cdfD = cdfU; family = missing, param = missing)

        if (ismissing(cdfU)); throw(ArgumentError("Cdf must be provided")); end
        if !(size(cdfU) == size(cdfD)); throw(ArgumentError("Upper and lower bounds must be the same size")); end

        return new(cdfU, cdfD, family, param)
    end
end

function (obj::copula)(X,Y, useInterp = false)             # Evaluating the copula

    #if (!issorted(X) || !issorted(Y)) throw(ArgumentError("cdf evaluation request must be given in assending order"));end

    nX, nY = size(obj.cdfU)

    if useInterp

        cdfU = interpolate(obj.cdfU, BSpline(Quadratic(Line(OnCell()))))
        cdfD = interpolate(obj.cdfD, BSpline(Quadratic(Line(OnCell()))))
        #println("interpolator used")
        x = X * (nX-1)  .+ 1;
        y = Y * (nY-1)  .+ 1;
        return interval(cdfD(x,y), cdfU(x,y));
    end

    xIndexLower = Int.(floor.(X .* (nX-1))) .+ 1
    yIndexLower = Int.(floor.(Y .* (nY-1))) .+ 1

    xIndexUpper = Int.(ceil.(X .* (nX-1)))  .+1
    yIndexUpper = Int.(ceil.(Y .* (nY-1)))  .+1

    #w = [opp(x,y) for x in X, y in Y]
    #m = [perf(x,y) for x in X, y in Y]

    #cdfLo = max.(obj.cdfD[xIndexLower, yIndexLower], w)
    #cdfHi = min.(obj.cdfU[xIndexUpper, yIndexUpper], m)

    cdfLo = obj.cdfD[xIndexLower, yIndexLower]
    cdfHi = obj.cdfU[xIndexUpper, yIndexUpper]

    return interval.(cdfLo, cdfHi)

end

cdf(C :: copula, X, Y) = C(X, Y)

function cdfU(C :: copula, X, Y)

    nX, nY = size(C.cdfU)

    xIndexLower = Int.(floor.(X .* (nX-1))) .+ 1
    yIndexLower = Int.(floor.(Y .* (nY-1))) .+ 1

    xIndexUpper = Int.(ceil.(X .* (nX-1)))  .+1
    yIndexUpper = Int.(ceil.(Y .* (nY-1)))  .+1

    #w = [opp(x,y) for x in X, y in Y]
    #m = [perf(x,y) for x in X, y in Y]
    
    #cdfLo = max.(C.cdfU[xIndexLower, yIndexLower], w)
    #cdfHi = min.(C.cdfU[xIndexUpper, yIndexUpper], m)

    cdfLo = C.cdfU[xIndexLower, yIndexLower]
    cdfHi = C.cdfU[xIndexUpper, yIndexUpper]

    return interval.(cdfLo, cdfHi)

end

function cdfD(C :: copula, X, Y)

    nX, nY = size(C.cdfD)

    xIndexLower = Int.(floor.(X .* (nX-1))) .+ 1
    yIndexLower = Int.(floor.(Y .* (nY-1))) .+ 1

    xIndexUpper = Int.(ceil.(X .* (nX-1)))  .+ 1
    yIndexUpper = Int.(ceil.(Y .* (nY-1)))  .+ 1
    
    #w = [opp(x,y) for x in X, y in Y]
    #m = [perf(x,y) for x in X, y in Y]
    
    #cdfLo = max.(C.cdfD[xIndexLower, yIndexLower], w)
    #cdfHi = min.(C.cdfD[xIndexUpper, yIndexUpper], m)

    cdfLo = C.cdfD[xIndexLower, yIndexLower]
    cdfHi = C.cdfD[xIndexUpper, yIndexUpper]

    return interval.(cdfLo, cdfHi)
end


function mass(C :: copula, x :: Interval, y :: Interval)

    C22hi = cdfU(C, x.hi, y.hi).hi
    C21hi = cdfD(C, x.hi, y.lo).lo
    C12hi = cdfD(C, x.lo, y.hi).lo
    C11hi = cdfU(C, x.lo, y.lo).hi

    C22lo = cdfD(C, x.hi, y.hi).lo
    C21lo = cdfU(C, x.hi, y.lo).hi
    C12lo = cdfU(C, x.lo, y.hi).hi
    C11lo = cdfD(C, x.lo, y.lo).lo

    Mhi = min(C22hi - C21hi - C12hi + C11hi, 1)
    Mlo = max(C22lo - C21lo - C12lo + C11lo, 0)

    return interval(Mlo, Mhi)

end

function massU(C :: copula, x :: Interval, y :: Interval)

    C22hi = cdfU(C, x.hi, y.hi).hi
    C21hi = cdfU(C, x.hi, y.lo).lo
    C12hi = cdfU(C, x.lo, y.hi).lo
    C11hi = cdfU(C, x.lo, y.lo).hi

    C22lo = cdfU(C, x.hi, y.hi).lo
    C21lo = cdfU(C, x.hi, y.lo).hi
    C12lo = cdfU(C, x.lo, y.hi).hi
    C11lo = cdfU(C, x.lo, y.lo).lo

    Mhi = min(C22hi - C21hi - C12hi + C11hi, 1)
    Mlo = max(C22lo - C21lo - C12lo + C11lo, 0)

    return interval(Mlo, Mhi)

end


function massD(C :: copula, x :: Interval, y :: Interval)

    C22hi = cdfD(C, x.hi, y.hi).hi
    C21hi = cdfD(C, x.hi, y.lo).lo
    C12hi = cdfD(C, x.lo, y.hi).lo
    C11hi = cdfD(C, x.lo, y.lo).hi

    C22lo = cdfD(C, x.hi, y.hi).lo
    C21lo = cdfD(C, x.hi, y.lo).hi
    C12lo = cdfD(C, x.lo, y.hi).hi
    C11lo = cdfD(C, x.lo, y.lo).lo

    Mhi = min(C22hi - C21hi - C12hi + C11hi, 1)
    Mlo = max(C22lo - C21lo - C12lo + C11lo, 0)

    return interval(Mlo, Mhi)

end

###
#   Copula generators for Archimedean (Frank and Clayton) copulas. Allows for easy and accurate copula generation 
#   at any dimension in terms of univariate functions (generator and inverse generator).
###

ClaGen(x, t = 1) = 1/t * (x ^(-t) -1 )       # Clayton copula generator
ClaInv(x, t = 1) = (1 + t * x) ^ (-1/t)      # Inverse generator

function FGen(x, s = 1)                      # Frank generator
    X1 = exp( -(x * s) ) -1
    X2 = exp( -s ) -1
    return -log( X1 / X2)
end

function FInv(x, s = 1)                      # Inverse
    X1 = exp( -x ) * (exp(-s) - 1)
    X2 = log(1 + X1)
    return - (X2 / s)
end


###
#   Copula functions and constructors
###
indep(x,y)          = x * y
perf(x,y)           = min(x,y)
opp(x,y)            = max(x+y -1,0)
Cla(x, y, t = 1)    = ClaInv( ClaGen(x, t) + ClaGen(y, t), t)
F(x,y,s = 1)        = FInv( FGen(x, s) + FGen(y, s), s)
function Gau(x,y,r=0)

    if x == 0; return 0; end
    if y == 0; return 0; end
    if x == 1; return y; end
    if y == 1; return x; end

    return bivariate_cdf(quantile.(Normal(),x),quantile.(Normal(),y), r)
end


###
#   Family constructors
###

function M()

    n = ProbabilityBoundsAnalysis.steps

    x = range(0, 1, length = n);
    cdf = [perf(xs,ys) for xs in x, ys in x]

    return copula(cdf, family = "M");
end

function W()

    n = ProbabilityBoundsAnalysis.steps

    x = range(0,1,length = n);
    cdf = [opp(xs,ys) for xs in x, ys in x]

    return copula(cdf, family = "W");
end

function πCop()
    
    n = ProbabilityBoundsAnalysis.steps

    x = range(0,1,length = n);
    cdf = [indep(xs, ys) for xs in x, ys in x]

    return copula(cdf, family = "π");
end
Pi() = πCop()

function Frank(s = 1)                       #   s is real; inf for perfect, 0 for indep, -inf for opposite

    if s == -Inf         # Limit should be set earlier
        C = W()
        return copula(C.cdfU, family = "Frank", param = "s = -Inf")
    end
    if s == 0
        C = πCop()
        return copula(C.cdfU,  family = "Frank", param =  "s = 1")
    end
    if s == Inf
        C = M()
        return copula(C.cdfU, family = "Frank",  param = "s = Inf")
    end

    n = ProbabilityBoundsAnalysis.steps
    x = range(0, 1, length = n);                     
    cdf = [F(xs, ys, s) for xs in x, ys in x];

    return copula(cdf, family = "Frank", param = "s = $s");
end

function Clayton(t = 0)                     #   t>-1; -1 for opposite, 0 for indep and inf for perfect

    if t == 0
        C = πCop()
        return copula(C.cdfU, family = "Clayton" , param = "t = 0")
    end
    if t == -1
        C = W()
        return copula(C.cdfU, family = "Clayton" , param = "t = -1")
    end
    if t == Inf                       # Limit should be set earlier
        C = M()
        return copula(C.cdfU, family = "Clayton" , param = "t = Ianf")
    end

    n = ProbabilityBoundsAnalysis.steps
    x = range(0,stop=1,length = n);
    cdf = [Cla(xs,ys, t) for xs in x, ys in x];

    return copula(cdf, family = "Clayton" , param = "t = $t");
end

function GauCopula(r = 0)     #   -1 <= r <=1 ; -1 for opposite, 1 for indep and 1 for perfect

    if r == 0; C = πCop(); return copula(C.cdfU, family = "Gaussian", param = "r = 0"); end
    if r == -1; C = W(); return copula(C.cdfU, family = "Gaussian", param = "r = -1"); end
    if r == 1; C = M(); return copula(C.cdfU, family = "Gaussian", param = "r = 1"); end

    n = ProbabilityBoundsAnalysis.steps
    x = range(0,stop = 1,length = n);

    cdf = [Gau(xs, ys ,r) for xs in x, ys in x];

    return copula(cdf, family = "Gaussian");
end

function Frechet()

    n = ProbabilityBoundsAnalysis.steps
    x = range(0,stop = 1,length = n);

    cdfU = [perf(xs,ys) for xs in x, ys in x]
    cdfD = [opp(xs,ys) for xs in x, ys in x]
    
    return copula(cdfU, cdfD, family = "Frechet");
end

function Base.show(io::IO, z::copula)

    statement1 = "Arb"
    statement2 = ""

    if !ismissing(z.family); statement1 = z.family;end
    if !ismissing(z.param); statement2 = z.param;end

    print(io, "Copula ~ $statement1($statement2)");
end

#=
function (obj::copula)(X :: Sampleable{Univariate},Y :: Sampleable{Univariate})             #Constructing a joint distribution using the copula
    return joint(X,Y,obj)
end
=#


#=
function piLoop() 
    x = range(0, 1, length = 200)
    [indep(xs,ys) for xs in x, ys in x]
end

indepVec(x,y) = x .* y

meshgrid(x,y) = (repeat(x',length(y),1),repeat(y,1,length(x)))

function piVec() 
    x = range(0, 1, length = 200)
    X, Y = meshgrid(x,x)
    cdf = indepVec(X[:], Y[:])
    cdf = reshape(cdf,200,200)
end
=#