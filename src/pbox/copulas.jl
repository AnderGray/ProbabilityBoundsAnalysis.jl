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

pyimport_conda("mpl_toolkits.mplot3d", "mpl_toolkits")
art3d = PyObject(PyPlot.art3D)
using3D()

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

    n = ProbabilityBoundsAnalysis.steps
    
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
        return copula(C.cdfU, family = "Frank", param = -Inf)
    end
    if s == 0
        C = πCop()
        return copula(C.cdfU,  family = "Frank", param =  1)
    end
    if s == Inf
        C = M()
        return copula(C.cdfU, family = "Frank",  param = Inf)
    end

    n = ProbabilityBoundsAnalysis.steps
    x = range(0, 1, length = n);                     
    cdf = [F(xs, ys, s) for xs in x, ys in x];

    return copula(cdf, family = "Frank", param = s);
end

function Frank(s :: Interval)

    if !(s ⊆ interval(-Inf,Inf)); throw(ArgumentError("coefficient must be ⊆ [-Inf, Inf]\n $r ⊄ [-Inf, Inf]"));end

    lb = Frank(s.lo);
    ub = Frank(s.hi);

    envelope = env(lb,ub);

    return copula(envelope.cdfU, envelope.cdfD, family= "Frank", param = s)
end


function Clayton(t = 0)                     #   t>-1; -1 for opposite, 0 for indep and inf for perfect

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

    n = ProbabilityBoundsAnalysis.steps
    x = range(0,stop=1,length = n);
    cdf = [Cla(xs,ys, t) for xs in x, ys in x];

    return copula(cdf, family = "Clayton" , param = t);
end

function Clayton(t :: Interval)

    if !(t ⊆ interval(-1,Inf)); throw(ArgumentError("coefficient must be ⊆ [-1, Inf]\n $r ⊄ [-1, Inf]"));end

    lb = Clayton(t.lo);
    ub = Clayton(t.hi);

    envelope = env(lb,ub);

    return copula(envelope.cdfU, envelope.cdfD, family= "Clayton", param = t)
end


function GauCopula(r = 0)     #   -1 <= r <=1 ; -1 for opposite, 1 for indep and 1 for perfect

    if r == 0; C = πCop(); return copula(C.cdfU, family = "Gaussian", param = 0); end
    if r == -1; C = W(); return copula(C.cdfU, family = "Gaussian", param =  -1); end
    if r == 1; C = M(); return copula(C.cdfU, family = "Gaussian", param = 1); end

    n = ProbabilityBoundsAnalysis.steps
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

    n = ProbabilityBoundsAnalysis.steps
    x = range(0,stop = 1,length = n);

    cdfU = [perf(xs,ys) for xs in x, ys in x]
    cdfD = [opp(xs,ys) for xs in x, ys in x]
    
    return copula(cdfU, cdfD, family = "Frechet");
end


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

    mass = 1/ProbabilityBoundsAnalysis.steps
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
    end

    print(io, "$statement0 ~ $statement1( $xShape(mean = $xMean, var = $xVar), $yShape(mean = $yMean, var = $yVar)$statement2)");
end

function plotStep(cop :: copula; name = missing, pn = 50, fontsize = 18, alpha = 0.3)

    AU = cop.cdfU
    AD = cop.cdfD
    m = size(AU)[1];
    if m < pn; ppn = m; else ppn = pn; end      #pn is the plot "sub-sample", we can't plot all elements

    x = y = range(0, stop = 1,length=m)
    #xgrid = repeat(x',ppn,1)
    #ygrid = repeat(y,1,ppn)

    nm = round(m/ppn); 

    x = x[1:Int(nm):end]
    oneIn = x[end] == 1
    
    if !oneIn; x = [x; 1]; end                             # We always want to include 1 in the plot

    zU = AU[1:Int(nm):end,1:Int(nm):end]
    zD = AD[1:Int(nm):end,1:Int(nm):end]

    if !oneIn;

        zU = [zU AU[end, 1:Int(nm):end]]
        zU = [zU' [AU[1:Int(nm):end, end]; 1]]

        zD = [zD AD[end, 1:Int(nm):end]]
        zD = [zD' [AD[1:Int(nm):end, end]; 1]]
    end

    if ismissing(name); fig = figure(figsize=(10,10)) else fig = figure(name,figsize=(10,10));end
    ax = fig.add_subplot(1,1,1,projection="3d")

    Ufaces = prepPolysU(zU,x);
    Dfaces = prepPolysD(zD,x);

    p3Ufaces = PyObject(art3d.Poly3DCollection(Ufaces, alpha=alpha))
    p3Dfaces = PyObject(art3d.Poly3DCollection(Dfaces, alpha=1))
    
    face_colorU = [1, 0, 0]
    face_colorD = [0, 0, 1]
    edge_color = [0.4, 0.4, 0.4]

    pycall(p3Ufaces.set_facecolor, PyAny, face_colorU)
    pycall(p3Ufaces.set_edgecolor, PyAny, edge_color)

    pycall(p3Dfaces.set_facecolor, PyAny, face_colorD)
    pycall(p3Dfaces.set_edgecolor, PyAny, edge_color)

    pycall(ax.add_collection3d, PyAny, p3Dfaces)
    pycall(ax.add_collection3d, PyAny, p3Ufaces)
    

end

function prepPolysU(us, x)

    ppn = length(x)

    p1 = [x[1]; x[1]; us[2, 2]]
    p2 = [x[1]; x[2]; us[2, 2]]
    p3 = [x[2]; x[2]; us[2, 2]]
    p4 = [x[2]; x[1]; us[2, 2]]
    tops = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]

    for i = 1:1:ppn-1
        for j = 1:1:ppn-1
            p1 = [x[i]; x[j]; us[i+1,j+1]]
            p2 = [x[i]; x[j+1]; us[i+1, j+1]]
            p3 = [x[i+1]; x[j+1]; us[i+1, j+1]]
            p4 = [x[i+1]; x[j]; us[i+1, j+1]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            tops = [tops; this]
        end
    end
    
    tops = tops[2:end]

    p1 = [x[1]; x[1]; us[1, 1]]
    p2 = [x[1]; x[1]; us[2, 2]]
    p3 = [x[1]; x[2]; us[2, 2]]
    p4 = [x[1]; x[2]; us[1, 1]]

    sides = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]

    for i = 1:1:ppn-1
        
        for j = i:1:ppn-1
            
            p1 = [x[i]; x[j]; us[i,j+1]]
            p2 = [x[i]; x[j]; us[i+1, j+1]]
            p3 = [x[i]; x[j+1]; us[i+1, j+1]]
            p4 = [x[i]; x[j+1]; us[i, j+1]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]
            
            p1 = [x[i]; x[j]; us[i+1,j]]
            p2 = [x[i]; x[j]; us[i+1, j+1]]
            p3 = [x[i+1]; x[j]; us[i+1, j+1]]
            p4 = [x[i+1]; x[j]; us[i+1, j]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]
            
            p1 = [x[j]; x[i]; us[j,i+1]]
            p2 = [x[j]; x[i]; us[j+1, i+1]]
            p3 = [x[j]; x[i+1]; us[j+1, i+1]]
            p4 = [x[j]; x[i+1]; us[j, i+1]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]
            
            p1 = [x[j]; x[i]; us[j+1,i]]
            p2 = [x[j]; x[i]; us[j+1, i+1]]
            p3 = [x[j+1]; x[i]; us[j+1, i+1]]
            p4 = [x[j+1]; x[i]; us[j+1, i]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]

        end
    end

    sides = unique(sides, dims=1)

    #return tops
    return [tops; sides]

end

function prepPolysD(us, x)

    ppn = length(x)

    p1 = [x[1]; x[1]; us[1, 1]]
    p2 = [x[1]; x[2]; us[1, 1]]
    p3 = [x[2]; x[2]; us[1, 1]]
    p4 = [x[2]; x[1]; us[1, 1]]
    tops = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]

    
    for i = 1:1:ppn-1
        for j = i:1:ppn-1
            p1 = [x[i]; x[j]; us[i+1, j]]
            p2 = [x[i]; x[j+1]; us[i+1, j]]
            p3 = [x[i+1]; x[j+1]; us[i+1, j]]
            p4 = [x[i+1]; x[j]; us[i+1, j]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            tops = [tops; this]
            
            p1 = [x[j]; x[i]; us[j, i+1]]
            p2 = [x[j]; x[i+1]; us[j, i+1]]
            p3 = [x[j+1]; x[i+1]; us[j, i+1]]
            p4 = [x[j+1]; x[i]; us[j, i+1]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            tops = [tops; this]

        end
    end
    
    p1 = [x[2]; x[2]; us[1, 1]]
    p2 = [x[2]; x[2]; us[2, 2]]
    p3 = [x[2]; x[3]; us[2, 2]]
    p4 = [x[2]; x[3]; us[1, 1]]

    sides = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
#=
    for i = 2:1:ppn-1
        for j = i:1:ppn-1
            
            p1 = [x[i-1]; x[j-1]; us[i-1,j]]
            p2 = [x[i-1]; x[j-1]; us[i, j]]
            p3 = [x[i-1]; x[j]; us[i, j]]
            p4 = [x[i-1]; x[j]; us[i-1, j]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]
            
            p1 = [x[i-1]; x[j-1]; us[i,j-1]]
            p2 = [x[i-1]; x[j-1]; us[i, j]]
            p3 = [x[i]; x[j-1]; us[i, j]]
            p4 = [x[i]; x[j-1]; us[i, j-1]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]
            
            p1 = [x[j-1]; x[i-1]; us[j-1,i]]
            p2 = [x[j-1]; x[i-1]; us[j, i]]
            p3 = [x[j-1]; x[i]; us[j, i]]
            p4 = [x[j-1]; x[i]; us[j-1, i]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]
            
            p1 = [x[j-1]; x[i-1]; us[j,i-1]]
            p2 = [x[j-1]; x[i-1]; us[j, i]]
            p3 = [x[j]; x[i-1]; us[j, i]]
            p4 = [x[j]; x[i-1]; us[j, i-1]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]

        end
    end
    =#

    
    for i = 1:1:ppn-1
        for j = i:1:ppn-1
            
            p1 = [x[i]; x[j]; us[i+1,j]]
            p2 = [x[i]; x[j]; us[i, j]]
            p3 = [x[i]; x[j+1]; us[i, j]]
            p4 = [x[i]; x[j+1]; us[i+1, j]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]
            #=
            p1 = [x[i+1]; x[j]; us[i,j+1]]
            p2 = [x[i+1]; x[j]; us[i+1, j]]
            p3 = [x[i]; x[j]; us[i+1, j]]
            p4 = [x[i]; x[j]; us[i, j+1]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]
            =#
            
            p1 = [x[j]; x[i]; us[j,i+1]]
            p2 = [x[j]; x[i]; us[j, i]]
            p3 = [x[j]; x[i+1]; us[j, i]]
            p4 = [x[j]; x[i+1]; us[j, i+1]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]
            
            #=
            p1 = [x[j]; x[i]; us[j+1,i]]
            p2 = [x[j]; x[i]; us[j+1, i+1]]
            p3 = [x[j+1]; x[i]; us[j+1, i+1]]
            p4 = [x[j+1]; x[i]; us[j+1, i]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]
            =#
        end
    end
    
    
    tops = unique(tops, dims=1)
    sides = unique(sides, dims=1)
    
    return [tops; sides]

end

function plot(x :: copula; name = missing, pn = 50, fontsize = 18, alpha = 0.7)

    AU = x.cdfU
    AD = x.cdfD
    m = size(AU)[1];
    if m < pn; ppn = m; else ppn = pn; end

    x = y = range(0, stop = 1,length=ppn)
    xgrid = repeat(x',ppn,1)
    ygrid = repeat(y,1,ppn)

    nm = round(m/ppn);

    zU = AU[1:Int(nm):end,1:Int(nm):end]
    zD = AD[1:Int(nm):end,1:Int(nm):end]

    if ismissing(name); fig = figure(figsize=(10,10)) else fig = figure(name,figsize=(10,10));end
    ax = fig.add_subplot(1,1,1,projection="3d")
    #ax = fig.add_subplot(2,1,1)
    
    plot_surface(xgrid, ygrid, zD, rstride=2,edgecolors="b", cstride=2, alpha=1, linewidth=1, cmap=ColorMap("Blues"))
    plot_surface(xgrid, ygrid, zU, rstride=2,edgecolors="r", cstride=2, alpha=alpha, linewidth=1, cmap=ColorMap("RdGy"))

    ax.view_init(45-27, 180+ 26)

    xlabel("y",fontsize = fontsize)
    ylabel("x",fontsize = fontsize)
    ax.zaxis.set_rotate_label(false);
    zlabel("C(x,y)", rotation = 0, fontsize = fontsize)
    #xticks(fontsize = fontsize); yticks(fontsize = fontsize);
    #PyPlot.title(title,fontsize=fontsize)
    tight_layout()
end



function plot(x :: bivpbox; name = missing, pn = 50, fontsize = 18, alpha = 0.7)

    AU = x.C.cdfU
    AD = x.C.cdfD
    m = size(AU)[1];
    if m < pn; ppn = m; else ppn = pn; end

    xus = x.marg1.u; xds = x.marg1.d;
    yus = x.marg2.u; yds = x.marg2.d;

    nm = round(m/ppn);

    xus = xus[1:Int(nm):end]; xds = xds[1:Int(nm):end]; 
    yus = yus[1:Int(nm):end]; yds = yds[1:Int(nm):end]; 

    zU = AU[1:Int(nm):end,1:Int(nm):end]
    zD = AD[1:Int(nm):end,1:Int(nm):end]

    uGridx, uGridy = (repeat(xus',length(yus),1),repeat(yus,1,length(xus)))
    dGridx, dGridy = (repeat(xds',length(yds),1),repeat(yds,1,length(xds)))

    if ismissing(name); fig = figure(figsize=(10,10)) else fig = figure(name,figsize=(10,10));end
    ax = fig.add_subplot(1,1,1,projection="3d")
    #ax = fig.add_subplot(2,1,1)
    
    plot_surface(dGridx, dGridy, zD, rstride=2,edgecolors="b", cstride=2, alpha=1, linewidth=1, cmap=ColorMap("Blues"))
    plot_surface(uGridx, uGridy, zU, rstride=2,edgecolors="r", cstride=2, alpha=alpha, linewidth=1, cmap=ColorMap("RdGy"))

    ax.view_init(45-27, 180+ 26)

    xlabel("y",fontsize = fontsize)
    ylabel("x",fontsize = fontsize)
    ax.zaxis.set_rotate_label(false);
    zlabel("H(x,y)", rotation = 0, fontsize = fontsize)
    #xticks(fontsize = fontsize); yticks(fontsize = fontsize);
    #PyPlot.title(title,fontsize=fontsize)
    tight_layout()
end




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