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

function (obj::copula)(X,Y, useInterp = false)             # Evaluating the copula (cdf)

    if X .< 0 || Y .< 0; return interval(0); end
    if X .> 1 || Y .> 1; return interval(1); end

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

(obj::copula)(X :: Array{T, 1}, Y :: Array{T, 1}, useInterp = false)   where T <: Real = [obj(i,j,useInterp) for i in X, j in Y]

(obj::copula)(X :: StepRangeLen, Y :: Array{T, 1}, useInterp = false) where T <: Real = [obj(i,j,useInterp) for i in X, j in Y]

(obj::copula)(X :: Array{T, 1}, Y :: StepRangeLen, useInterp = false)  where T <: Real = [obj(i,j,useInterp) for i in X, j in Y]

(obj::copula)(X :: StepRangeLen, Y :: StepRangeLen, useInterp = false) = [obj(i,j,useInterp) for i in X, j in Y]

(obj::copula)(X :: Real, Y :: Array{T, 1}, useInterp = false)   where T <: Real = [obj(X,j,useInterp) for j in Y]

(obj::copula)(X :: Array{T, 1}, Y :: Real, useInterp = false)   where T <: Real = [obj(i,Y,useInterp) for i in X]

(obj::copula)(X :: Real, Y :: StepRangeLen, useInterp = false)  = [obj(X,j,useInterp) for j in Y]

(obj::copula)(X :: StepRangeLen, Y :: Real, useInterp = false)  = [obj(i,Y,useInterp) for i in X]

function (obj::copula)(X :: Interval, Y::Interval, useInterp = false)

    cdfLo = obj(X.lo,Y.lo, useInterp).lo
    cdfHi = obj(X.hi,Y.hi, useInterp).hi

    return interval(cdfLo, cdfHi)

end

function (obj::copula)(X :: Interval, Y::Real, useInterp = false)
    Y = interval(Y);
    return obj(X,Y,useInterp)
end

function (obj::copula)(X :: Real, Y::Interval, useInterp = false)
    X = interval(X);
    return obj(X,Y,useInterp)
end

function (obj::copula)(X :: pbox,Y :: pbox)             #Constructing a joint distribution using the copula
    return bivpbox(X, Y, obj)
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

function conditionalX(C :: AbstractCopula, xVal :: Real)

    Nsx, Nsy = size(C.cdfU)

    n = parametersPBA.steps

    yGridU = ProbabilityBoundsAnalysis.ii();
    yGridD = ProbabilityBoundsAnalysis.jj();

    diff = 1/Nsy * 2

    zU = (left.(C(xVal + diff/2, yGridU)  - left.(C(xVal - diff/2, yGridU))))/diff
    zD = (right.(C(xVal + diff/2, yGridD) - right.(C(xVal - diff/2, yGridD))))/diff

    ###
    #   zU are the cdf values, and yGridU/D are the physical values, which is a uniform grid.
    #   We need to "invert" it so that zU is a uniform grid.
    ###

    inverseBox = pbox(zU,zD)

    us = left.(cdf.(inverseBox,yGridU))
    ds = right.(cdf.(inverseBox,yGridD))

    return pbox(us,ds)

end

function conditionalY(C :: AbstractCopula, yVal :: Real)

    Nsx, Nsy = size(C.cdfU)
    xGridU = ii();
    xGridD = jj();

    diff = 1/Nsx * 2

    zU = (left.(C(xGridU, yVal + diff/2) - left.(C(xGridU, yVal - diff/2))))/diff
    zD = (right.(C(xGridD,yVal + diff/2) - right.(C(xGridD,yVal - diff/2))))/diff

    ###
    #   zU are the cdf values, and yGridU/D are the physical values, which is a uniform grid.
    #   We need to "invert" it so that zU is a uniform grid.
    ###


    inverseBox = pbox(zU,zD)

    us = left.(cdf.(inverseBox,xGridU))
    ds = right.(cdf.(inverseBox,xGridD))

    return pbox(us,ds)


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

# Copula functions

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
# Copula bounds from Kendal τ and Spearman ρ
# BOUNDS ON BIVARIATE DISTRIBUTION FUNCTIONS WITH GIVEN MARGINS AND MEASURES OF ASSOCIATION, Roger B. Nelsen et al.
###

#kenLB(x, y, τ) = max(0, x+y-1, 0.5*( (x+y)-sqrt( (x-y)^2)+1-τ ) )
kenUB(x, y, τ) = min(x, y, 0.5*( (x+y-1) + sqrt( (x+y-1)^2 +1+τ )) )
kenLB(x, y, τ) = x - kenUB(x, 1 - y, - τ)

#spearϕ(a, b) = 1/6 * ( (max(0, 9*b + 3*sqrt( 9*b^2 - 3*a^6)))^(1/3) + ( max(0,9*b - 3*sqrt(9*b^2 - 3*a^6)))^(1/3) )

function spearϕ(a, b)
    A = 9*b
    B = max(9*b^2 - 3*a^6, 0)
    C = (max(0, A + 3*sqrt(B)))^(1/3)
    D = (max(0, A - 3*sqrt(B)))^(1/3)
    return 1/6 * (C + D)
end

#spearLB(x, y, ρ) = max(0, x + y - 1, (x + y)/2 - spearϕ(x - y, 1 - ρ))
spearUB(x, y, ρ) = max( opp(x,y), min(x, y, (x + y -1)/2 + spearϕ(x + y - 1, 1 + ρ)))
spearLB(x, y, τ) = x - spearUB(x, 1 - y, - τ)

###
#   Family constructors
###

function M()

    n = parametersPBA.steps

    x = range(0, 1, length = n);
    cdf = [perf(xs,ys) for xs in x, ys in x]

    return copula(cdf, family = "M");
end

function W()

    n = parametersPBA.steps

    x = range(0,1,length = n);
    cdf = [opp(xs,ys) for xs in x, ys in x]

    return copula(cdf, family = "W");
end

function πCop()

    n = parametersPBA.steps

    x = range(0,1,length = n);
    cdf = [indep(xs, ys) for xs in x, ys in x]

    return copula(cdf, family = "π");
end
Pi() = πCop()

function Frank(s = 1)                       #   s is real; inf for perfect, 0 for indep, -inf for opposite

    if !(s ∈ interval(-Inf,Inf)); throw(ArgumentError("coefficient must be ∈ [-Inf, Inf]\n $s ∉ [-Inf, Inf]"));end

    if s == -Inf         # Limit should be set earlier
        C = W()
        return copula(C.cdfU, family = "Frank", param = -Inf)
    end
    if s == 0
        C = πCop()
        return copula(C.cdfU,  family = "Frank", param =  0)
    end
    if s == Inf
        C = M()
        return copula(C.cdfU, family = "Frank",  param = Inf)
    end

    n = parametersPBA.steps
    x = range(0, 1, length = n);
    cdf = [F(xs, ys, s) for xs in x, ys in x];

    return copula(cdf, family = "Frank", param = s);
end

function Frank(s :: Interval)

    if !(s ⊆ interval(-Inf,Inf)); throw(ArgumentError("coefficient must be ⊆ [-Inf, Inf]\n $s ⊄ [-Inf, Inf]"));end

    lb = Frank(s.lo);
    ub = Frank(s.hi);

    envelope = env(lb,ub);

    return copula(envelope.cdfU, envelope.cdfD, family= "Frank", param = s)
end


function Clayton(t = 0)                     #   t>-1; -1 for opposite, 0 for indep and inf for perfect

    if !(t ∈ interval(-1,Inf)); throw(ArgumentError("coefficient must be ∈ [-1, Inf]\n $t ∉ [-1, Inf]"));end

    if t == 0
        C = πCop()
        return copula(C.cdfU, family = "Clayton" , param = 0)
    end
    if t == -1
        C = W()
        return copula(C.cdfU, family = "Clayton" , param = -1)
    end
    if t == Inf                       # Limit should be set earlier
        C = M()
        return copula(C.cdfU, family = "Clayton" , param = Inf)
    end

    n = parametersPBA.steps
    x = range(0,stop=1,length = n);
    cdf = [Cla(xs,ys, t) for xs in x, ys in x];

    return copula(cdf, family = "Clayton" , param = t);
end

function Clayton(t :: Interval)

    if !(t ⊆ interval(-1,Inf)); throw(ArgumentError("coefficient must be ⊆ [-1, Inf]\n $t ⊄ [-1, Inf]"));end

    lb = Clayton(t.lo);
    ub = Clayton(t.hi);

    envelope = env(lb,ub);

    return copula(envelope.cdfU, envelope.cdfD, family= "Clayton", param = t)
end


function GauCopula(r = 0)     #   -1 <= r <=1 ; -1 for opposite, 1 for indep and 1 for perfect

    if !(r ∈ interval(-1,1)); throw(ArgumentError("correlation coefficient must be ∈ [-1, 1]\n $r ∉ [-1, 1]"));end

    if r == 0; C = πCop(); return copula(C.cdfU, family = "Gaussian", param = 0); end
    if r == -1; C = W(); return copula(C.cdfU, family = "Gaussian", param =  -1); end
    if r == 1; C = M(); return copula(C.cdfU, family = "Gaussian", param = 1); end

    n = parametersPBA.steps
    x = range(0,stop = 1,length = n);

    cdf = [Gau(xs, ys ,r) for xs in x, ys in x];

    return copula(cdf, family = "Gaussian", param = r);
end

function GauCopula(r :: Interval)

    if !(r ⊆ interval(-1,1)); throw(ArgumentError("correlation coefficient must be ⊆ [-1, 1]\n $r ⊄ [-1, 1]"));end

    lb = GauCopula(r.lo);
    ub = GauCopula(r.hi);

    envelope = env(lb,ub);

    return copula(envelope.cdfU, envelope.cdfD, family= "Gaussian", param = r)
end

function Frechet()

    n = parametersPBA.steps
    x = range(0,stop = 1,length = n);

    cdfU = [perf(xs,ys) for xs in x, ys in x]
    cdfD = [opp(xs,ys) for xs in x, ys in x]

    return copula(cdfU, cdfD, family = "Frechet");
end


function τCopula( τ = 0 )   # Imprecise copula from kendal tau

    n = parametersPBA.steps
    x = range(0,stop = 1,length = n);

    cdfU = [kenUB(xs,ys,τ) for xs in x, ys in x]
    cdfD = [kenLB(xs,ys,τ) for xs in x, ys in x]

    return copula(cdfU, cdfD, family = "Kendal", param = τ)

end

KendalCopula(τ = 0) = τCopula(τ)

function ρCopula( ρ = 0 ) # Imprecise copula from Spearman rho

    n = parametersPBA.steps
    x = range(0,stop = 1,length = n);

    cdfU = [spearUB(xs,ys,ρ) for xs in x, ys in x]
    cdfD = [spearLB(xs,ys,ρ) for xs in x, ys in x]

    return copula(cdfU, cdfD, family = "Spearman", param = ρ)
end

SpearmanCopula(ρ = 0 ) = ρCopula( ρ )

function env(x::copula, y::copula)

    sameU = size(x.cdfU) == size(y.cdfU);
    sameD = size(x.cdfU) == size(y.cdfU);

    if !sameU || !sameD; throw(ArgumentError("descritization of provided copulas differ")); end

    newU = max.(x.cdfU, y.cdfU);
    newD = min.(x.cdfD, y.cdfD);

    if x.family == y.family; newFam = x.family; else newFam = missing; end
    return copula(newU, newD, family = newFam)
end

###
#   Copula sampling
###
function sampleCond(C :: AbstractCopula, N = 1)

    if !(C.cdfU == C.cdfU); throw(ArgumentError("Can only sample from precise copula")); end

    useInterp = false;      # Set true to use interpolator
    if ismissing(C.family); family = 0; useInterp = true else family = C.family end

    if (family == "Gaussian") return CholeskyGaussian(N, Float64(C.param)) ;end  # Use Cholesky decompostition of the Cov matrix for gaussian copula sampling

    n = size(C.cdfD)[1]

    x = rand(N);    y = rand(N);
    ux = x;         uy = zeros(N);
    m = n;

    if family == "M"; return hcat(ux,ux); end
    if family == "W"; return hcat(ux,1 .- ux); end
    if family == "π"; return hcat(ux,y); end

    uy = mid.([cut(conditionalY(C, x[i]), y[i]) for i =1:N])

    return hcat(ux,uy);
end

function sample(C :: AbstractCopula, N = 1)

    if !(C.cdfU == C.cdfU); throw(ArgumentError("Can only sample from precise copula")); end

    useInterp = false;      # Set true to use interpolator
    if ismissing(C.family); family = 0; useInterp = true else family = C.family end

    if (family == "Gaussian") return CholeskyGaussian(N, Float64(C.param)) ;end  # Use Cholesky decompostition of the Cov matrix for gaussian copula sampling

    n,n2 = size(C.cdfD)


    x = rand(N);    y = rand(N);
    ux = x;         uy = zeros(N);
    m = n;

    if family == "M"; return hcat(ux,ux); end
    if family == "W"; return hcat(ux,1 .- ux); end
    if family == "π"; return hcat(ux,y); end

    #uy = mid.([cut(conditionalY(C, x[i]), y[i]) for i =1:N])

    uMasses = C.cdfU[2:end, 2:end] + C.cdfU[1:end-1, 1:end-1] - C.cdfU[1:end-1, 2:end] - C.cdfU[2:end, 1:end-1]
    dMasses = C.cdfD[2:end, 2:end] + C.cdfD[1:end-1, 1:end-1] - C.cdfD[1:end-1, 2:end] - C.cdfD[2:end, 1:end-1]

    is = range(0,1,length = n); js = range(0,1,length = n2)
    #boxesU = interval.(is[1:end-1], is[2:end]) .× interval.(js[1:end-1], js[2:end])
    boxesU = [interval(is[i], is[i+1]) × interval(js[j], js[j+1]) for i = 1:n-1, j = 1:n2-1]

    uMass = uMasses[:]
    boxesU = boxesU[:]

    if !(sum(uMass) == 1); uMass = uMass./sum(uMass);end
    cdfUs = cumsum(uMass)
    ux = interval.(zeros(N), zeros(N)); uy = interval.(zeros(N), zeros(N))
    for i =1:N
        index = findfirst(x[i] .< cdfUs)
        ux[i] = boxesU[index][1];
        uy[i] = boxesU[index][2];
    end

    return hcat(ux,uy);
end



function CholeskyGaussian(N = 1, correlation = 0)

    Cov = ones(2,2); Cov[1,2] = Cov[2,1] = correlation;
    a = cholesky(Cov).L;
    z = transpose(rand(Normal(),N,2))
    x = a * z;
    u = transpose(cdf.(Normal(),x))

    return hcat(u[:,1],u[:,2])
end


###
#   Copula rotations. Rotates the mass by 90° (makes M -> W), 180° or 270°
###
#=
# C⁹⁰(u1, u2) = u1 - C(1-u2, u1)
function rotate90(x :: copula)
    cdfU = x.cdfU; cdfD = x.cdfD

    n1,n2 = size(cdfU);

    u1 = range(0, 1, length = n1);
    u2 = range(0, 1, length = n2);

    cdfU = u1 .- right.(cdfU(x, 1 .- u2, u1))

end

# C180(u1,u2) = u1 + u2 + C(1 - u1, 1 - u2) -1
function rotate180(x :: copula)
    cdfU = x.cdfU; cdfD = x.cdfD

    n1,n2 = size(cdfU);

    u1 = range(0, 1, length = n1);

end

# C270(u1,u2) = u2 - C(u2, 1 - u1)
function rotate270(x :: copula)
    cdfU = x.cdfU; cdfD = x.cdfD

    n1,n2 = size(cdfU);

    u1 = range(0, 1, length = n1);

end
=#
###
#   Bivaraite pbox
###
mutable struct bivpbox <: AbstractJoint

    marg1 :: pbox
    marg2 :: pbox
    C  :: AbstractCopula

    function bivpbox(marginal1, marginal2, copula = πCop())

        return new(marginal1, marginal2, copula)

    end
end

function (obj::bivpbox)(X :: Real, Y :: Real, useInterp = false)             # Evaluating the cdf, sklar:  C(F(x),G(y))

    cdfsX = cdf(obj.marg1, X)
    cdfsY = cdf(obj.marg2, Y)

    return obj.C(cdfsX, cdfsY, useInterp)

end

function (obj::bivpbox)(X, Y, useInterp = false)             # Evaluating the cdf, sklar:  C(F(x),G(y))

    cdfsX = [cdf(obj.marg1, x) for x in X]
    cdfsY = [cdf(obj.marg2, y) for y in Y]

    return [obj.C(is, js, useInterp) for is in cdfsX, js in cdfsY]

end

cdf(Biv :: bivpbox, X, Y, useInterp = false) = Biv(X,Y,useInterp)

function mass(Biv :: bivpbox, X :: Interval, Y :: Interval)


    H11 = Biv(X.lo,Y.lo)
    H12 = Biv(X.lo,Y.hi)
    H21 = Biv(X.hi,Y.lo)
    H22 = Biv(X.hi,Y.hi)

    Mhi = min(H22.hi - H21.lo - H12.lo + H11.hi, 1)
    Mlo = max(H22.lo - H21.hi - H12.hi + H11.lo, 0)

    return interval(Mlo, Mhi)
end

function conditionalX(J :: bivpbox, xVal :: Real)

    #Nsx, Nsy = size(J.C.cdfU)
    yGridU = J.marg2.u;
    yGridD = J.marg2.d;

    diff = maximum(yGridU[2:end] .- yGridU[1:end-1])

    zU = (left.(J(xVal + diff/2, yGridU) - left.(J(xVal - diff/2, yGridU))))/diff
    zD = (right.(J(xVal + diff/2, yGridD) - right.(J(xVal - diff/2, yGridD))))/diff

    uIndex = findfirst(xVal .< yGridU);
    dIndex = findfirst(xVal .< yGridD);

    mass = 1/parametersPBA.steps
    densityU = mass/(yGridU[uIndex] - yGridU[uIndex-1])
    densityD = mass/(yGridD[dIndex] - yGridD[dIndex-1])

    zU = zU ./densityU;     zD = zD ./densityD;

    inverseBox = pbox(zU,zD)

    us = left.(cdf.(inverseBox,yGridU))
    ds = right.(cdf.(inverseBox,yGridD))


    return pbox(us,ds)

end


function conditionalY(J :: bivpbox, yVal :: Real)

    Nsx, Nsy = size(J.C.cdfU)

    xGridU = J.marg2.u;
    xGridD = J.marg2.d;

    diff = maximum(xGridU[2:end] .- xGridU[1:end-1])
    zU = (left.(J(xGridU, yVal + diff/2) - left.(J(xGridU, yVal - diff/2))))/diff
    zD = (right.(J(xGridD,yVal + diff/2) - right.(J(xGridD,yVal - diff/2))))/diff

    return pbox(zU, zD)

end

###
#   Sampling biv pbox
###
function sampleCond(J :: bivpbox, N = 1)

    C = J.C; x = J.marg1; y = J.marg2
    copSamples = sampleCond(C,N)

    ux = cut.(x, copSamples[:,1])
    uy = cut.(y, copSamples[:,2])

    return hcat(ux,uy);
end

function sample(J :: bivpbox, N = 1)

    C = J.C; x = J.marg1; y = J.marg2
    copSamples = sample(C,N)

    ux = cut.(x, copSamples[:,1])
    uy = cut.(y, copSamples[:,2])

    return hcat(ux,uy);
end



###
#   Plotting and show functions
###

function Base.show(io::IO, z::copula)

    statement0 = "Copula"
    statement1 = "Arb"
    statement2 = ""

    if !(z.cdfU == z.cdfD); statement0 = "Imp Copula";end

    if !ismissing(z.family); statement1 = z.family;end
    if !ismissing(z.param)
        if z.family == "Gaussian"; statement2 = "r = $(z.param)"; end
        if z.family == "Frank"; statement2 = "s = $(z.param)"; end
        if z.family == "Clayton"; statement2 = "t = $(z.param)"; end
    end

    print(io, "$statement0 ~ $statement1($statement2)");
end


function Base.show(io::IO, z::bivpbox)

    statement0 = "BivPbox"
    statement1 = "Arb"
    statement2 = ""

    range = "[$(left(z.marg1)), $(right(z.marg1))] x [$(left(z.marg2)), $(right(z.marg2))]"

    xShape = z.marg1.shape; yShape = z.marg2.shape;

    xMean = mean(z.marg1); yMean = mean(z.marg2);
    xVar = var(z.marg1); yVar = var(z.marg2);

    if xMean.lo == xMean.hi; xMean = xMean.hi; end
    if yMean.lo == yMean.hi; yMean = yMean.hi; end
    if xVar.lo == xVar.hi; xVar = xVar.hi; end
    if yVar.lo == yVar.hi; yVar = yVar.hi; end


    if !ismissing(z.C.family); statement1 = z.C.family;end
    if !ismissing(z.C.param)
        if z.C.family == "Gaussian"; statement2 = "; r = $(z.C.param)"; end
        if z.C.family == "Frank"; statement2 = "; s = $(z.C.param)"; end
        if z.C.family == "Clayton"; statement2 = "; t = $(z.C.param)"; end
        if z.C.family == "Spearman"; statement2 = "; ρ = $(z.C.param)"; end
        if z.C.family == "Kendal"; statement2 = "; τ = $(z.C.param)"; end
    end

    print(io, "$statement0 ~ $statement1( $xShape(mean = $xMean, var = $xVar), $yShape(mean = $yMean, var = $yVar)$statement2)");
end
