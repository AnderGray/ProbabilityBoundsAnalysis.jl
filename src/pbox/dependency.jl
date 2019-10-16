######
# This file is part of the pba.jl package.
#
#   Definition of bivariate copula class plus methods
#
#           University of Liverpool
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
###

###
#
#   Known bugs:
#                   ->  The differentiation for the sampling should not depend on n. May request values x>1 or x<0
#                   ->  eps() should be used for differentiation. Liams code.
#                   ->  Linear interpolator not working...
#                   ->  Bilinear interpolartor? 
#
###

using Distributions
using PyPlot
using3D()
include("NormalDistribution.jl")

global  n = 200;        # Number of descretizations
global  pn = 200;        # Number of descretization for plotting

abstract type AbstractCopula <: Real end
abstract type AbstractJoint <: Real end

mutable struct copula <: AbstractCopula

    cdf ::  Array{Float64,2}
    density :: Array{Float64,2}
    func :: Union{Function,Missing}
    param :: Union{Float64,Missing}

    function copula(cdf = missing; density = missing, func = missing, param = missing)

        if (ismissing(cdf) && ismissing(density)) throw(ArgumentError("Cdf and/or density must be provided")) end
        if (ismissing(cdf)); cdf = calcCdf(density); end
        if (ismissing(density)); density = calcDensity(cdf);end
        new(cdf,density, func, param)

    end
end

function (obj::copula)(X,Y)             # Evaluating the copula

    if (!issorted(X) || !issorted(Y)) throw(ArgumentError("cdf evaluation request must be given in assending order"));end

    if !ismissing(obj.func)
        if !ismissing(obj.param)
            return obj.func(X,Y,obj.param)
        end
        return obj.func(X,Y)
    end

    Indx = Int.(round.(collect(X*n)));   Indy = Int.(round.(collect(Y*n)));

    if (Indx[1] == 0) if (Indx[2] == 1)   deleteat!(Indx,1); else Indx[1] = 1; end; end
    if (Indy[1] == 0) if (Indy[2] == 1)   deleteat!(Indy,1); else Indy[1] = 1; end; end
    return obj.cdf[Indx,Indy];
end

function sample(C :: AbstractCopula, N = 1)
    x = rand(N);    y = rand(N);
    ux = x;         uy = zeros(N);
    m = n;

    if C.func   == perf return hcat(ux,ux); end
    if C.func   == opp return hcat(ux,1 .- ux); end
    if C.func   == indep return hcat(ux,y); end

    # if (!ismissing(C.func) && n < 10^5) m = 10^5; end  # If function is provided more acurate sampling, otherwise must use interpolator
    e = 1/m;        ygrid = range(0,1,length = m);
    #if !ismissing(C.func) e = eps();end
    for i=1:N
        conditional = (C(x[i] + e,ygrid) - C(x[i],ygrid))./e
        uy[i] = findlast(conditional .<= y[i]) / m;
    end
    return hcat(ux,uy);
end

function conditional(C :: AbstractCopula, xVal :: Real; plot = false)

    e = 1/n;        ygrid = range(0,1,length = n);
    if !ismissing(C.func) e = eps(); end

    conditional = (C(xVal + e/2,ygrid) - C(xVal - e/2,ygrid))./e;

    if (plot) plot2(C.cdf, conditional, xVal) end

    return conditional
end

function calcCdf(density :: Array{<:Real,2})
    # Integral of the density       Liams code
    return density
end

function calcDensity(cdf :: Array{<:Real,2})
    # Derriviative of the cdf (Derrivate of the Derrivate)
    return cdf
end

function M()
    x = y = range(0,1,length = n);
    return copula(perf(x,y), func = perf);
end

function W()
    x = y = range(0,1,length = n);
    return copula(opp(x,y), func = opp);
end

function π()
    x = y = range(0,1,length = n);
    return copula(indep(x,y), func = indep);
end
Pi() = π()

function Frank(s = 200)                                     #   Frank copula
    x = y = range(0,1,length = n);                          #   s>0; 0 for perfect, 1 for indep, inf for oposite
    return copula(F(x,y,s), func = F, param = s);
end

function Clayton(t = -0.5)                                  #   Clayton copula
    x = y = range(0,1,length = n);                          #   t>-1; -1 for opposite, 0 for indep and inf for perfect
    return copula(Cla(x,y,t), func = Cla, param = t);
end

function Gaussian(corr = 0)
    x = y = range(0,1,length = n);
    return copula(Gau(x,y,corr), func = Gau, param = corr);
end

indep(X, Y) = [x*y for x in X, y in Y];
perf(X, Y)  = [min(x,y) for x in X, y in Y];
opp(X, Y)   = [max(x+y-1,0) for x in X, y in Y];
F(X,Y,s = 200)      = [log(1+(s^x-1)*(s^y-1)/(s-1))/log(s) for x in X, y in Y]
Cla(X,Y, t = - 0.5) = [max((x^(-t)+y^(-t)-1)^(-1/t),0) for x in X, y in Y]
Gau(X,Y, corr = 0)  = [bivariate_cdf(quantile.(Normal(),x),quantile.(Normal(),y), corr) for x in X, y in Y];

mvNormCdf(X,Y,cor) = [bivariate_cdf(x,y,cor) for x in X, y in Y];

wait_for_key() = (print(stdout, "press enter to continue"); read(stdin, 1); nothing);

function samplePlot(C :: AbstractCopula, N = 1)

    x = rand(N);    y = rand(N);
    ux = x;         uy = zeros(N);
    m = n;
    e = 1/m;        ygrid = range(0,1,length = m);
    for i=1:N
        conditional = (C(x[i] + e/2,ygrid) - C(x[i] - e/2,ygrid))./e
        uy[i] = findlast(conditional .<= y[i]) / m;
        plot2(C.cdf,conditional,ux[i], uy[i], y[i])
        println(hcat(ux[i],uy[i]))
        wait_for_key()
    end
    return hcat(ux,uy);
end

plotDensity(x :: AbstractCopula) = plot1(x.density);
plotCdf(x::AbstractCopula) = plot1(x.cdf);

function plot1(A :: Array{<:Real,2})

#   https://github.com/gizmaa/Julia_Examples/blob/master/pyplot_surfaceplot.jl
    x = y = range(0,1,length=pn)
    xgrid = repeat(x',pn,1)
    ygrid = repeat(y,1,pn)
    z = A[1:Int(n/pn):end,1:Int(n/pn):end]

    fig = figure("SurfacPlots",figsize=(10,10))
    ax = fig.add_subplot(2,1,1,projection="3d")
    plot_surface(xgrid, ygrid, z, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("gray"), alpha=0.8, linewidth=0.25)
    xlabel("X")
    ylabel("Y")
    PyPlot.title("Surface Plot")

    subplot(212)
    ax = fig.add_subplot(2,1,2)
    cp = contour(xgrid, ygrid, z, colors="black")
    ax.clabel(cp, inline=1, fontsize=10)
    xlabel("X")
    ylabel("Y")
    PyPlot.title("Contour Plot")
    tight_layout()

end

function plot2(A :: Array{<:Real,2}, cond:: Array{<:Real,1}, ux = missing, uy = missing, uuy = missing)

    PyPlot.clf()
    x = y = range(0,1,length=pn)
    xgrid = repeat(x',pn,1)
    ygrid = repeat(y,1,pn)
    z = A[1:Int(n/pn):end, 1:Int(n/pn):end]
    conditional = cond[1:Int(n/pn):end];

    fig = figure("ConditionalPlots",figsize=(10,10))
    ax = fig.add_subplot(2,1,1)
    cp = contour(xgrid, ygrid, z)
    ax.clabel(cp, inline=1, fontsize=10)
    xlabel("X")
    ylabel("Y")

    if !ismissing(ux)
        xx = ones(pn) * ux
        plot(xx,x, color = :red,linewidth=2,label = "Rand1: $ux");
    end

    subplot(212)
    ax = fig.add_subplot(2,1,2)
    plot(y, conditional, color = :black)
    if !ismissing(uy) && !ismissing(uuy)
        nums1 = Int(round(pn*uy))
        nums2 = Int(round(pn*uuy))
        yy = ones(nums1) * uuy
        uyy = ones(nums2) * uy
        plot(y[1:nums1],yy,color = :red,linewidth=2, label= "Rand2: $uy");
        plot(uyy,y[1:nums2],color = :blue,linewidth=2);
    end
    xlabel("Y")
    ylabel("conditional CDF")
    tight_layout()
end

#=
function linInterp(A :: Array{<:Real,1}, x :: Real)

    #uy[i] = findlast(conditional .<= y[i]) / m;

    Ind = findlast(A .< x);
    numel = length(A);
    x0 = A[Ind];    x1 = A[Ind+1];
    y0 = Ind/numel; y1 = (Ind+1)/numel;

    #return (y0 * (x1-x) + y1 * (x - x0))/(x1-x0)
    return y0 * ((x1-x)/(x1-x0)) + y1 * ((x-x0)/(x1-x0))
    #sepX = 1/length(A);   # Uniform Grid and A[1] == 0 assumed

end

=#
#=
mutable struct joint <: AbstractCopula

    marginal1
    marginal2
    copula
    # The above properties will be enough...

    cdf ::  Array{Float64,2}
    density :: Array{Float64,2}
    func :: Union{Function,Missing}
    param :: Union{Float64,Missing}

    function joint(cdf = missing; density = missing, func = missing, param = missing)

        if (ismissing(cdf) && ismissing(density)) throw(ArgumentError("Cdf and/or density must be provided")) end
        if (ismissing(cdf)); cdf = calcCdf(density); end
        if (ismissing(density)); density = calcDensity(cdf);end
        new(cdf,density, func, param)

    end
end

function bilinearInterp(A :: Array{<:Real,2}, x:: Array{<:Real,1})

    if (length(x) != ndims(A)) throw(ArgumentError("Dimension of function and point must both be 2")); end
    if (length(x) != 2)        throw(ArgumentError("Dimension point must both be 2")); end
    if (ndims(A) != 2)        throw(ArgumentError("Dimension function must both be 2")); end

end

=#
