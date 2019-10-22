######
# This file is part of the pba.jl package.
#
#   Definition of bivariate copula class, joint distribution class and methods
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
#
#   For algebra of random variables, including algebra of covaraince See https://en.wikipedia.org/wiki/Algebra_of_random_variables
#       -> The above link is also useful for moment propagation. Could you perform covariance algebra when the covarainces are partially known
#
#
#   One model for dependency trancking, is to assume a copula family, and compute the covariance of the results of opartaions with their inputs
#   using the arithemtic of covaraince described above. For this we would need an interval covaraince algebra, and a way of performing convolutions with it.
#
###

###
#
#   Known bugs:
#                   ->  The differentiation for the sampling should not depend on n. May request values x>1 or x<0  (Solved)
#                   ->  eps() should be used for differentiation. Liams code. (Solved)
#                   ->  Linear interpolator not working... (Solved)
#                   ->  Bilinear interpolartor? (Solved)
#
###

using Distributions
using PyPlot
using PyCall
using LinearAlgebra
using Interpolations
using3D()
include("NormalDistribution.jl")

global  n  = 200;        # Number of descretizations
global  pn = 200;        # Number of descretization for plotting
global bOt = 0.001;              # smallest quamtile to use if left tail is unbounded
global tOp = 0.999;              # largest quamtile to use if right tail is unbounded

abstract type AbstractCopula <: Real end
abstract type AbstractJoint <: Real end

mutable struct copula <: AbstractCopula

    cdf ::  Array{Float64,2}
    density :: Union{Array{Float64,2},Missing}
    func :: Union{Function,Missing}
    param :: Union{Float64,Missing}

    function copula(cdf = missing; density = missing, func = missing, param = missing)

        if (ismissing(cdf) && ismissing(density)) throw(ArgumentError("Cdf and/or density must be provided")) end
        if (ismissing(cdf)); cdf = calcCdf(density); end
        C = new(cdf,density,func,param)
        if (ismissing(density)); density = calcDensity(C); end

        return new(cdf,density, func, param)

    end
end

function (obj::copula)(X,Y, useInterp = false)             # Evaluating the copula

    #if (!issorted(X) || !issorted(Y)) throw(ArgumentError("cdf evaluation request must be given in assending order"));end

    if !useInterp
        if !ismissing(obj.func)
            if !ismissing(obj.param)
                return obj.func(X,Y,obj.param)
            end
            return obj.func(X,Y)
        end
    end

    cdf = interpolate(obj.cdf, BSpline(Quadratic(Line(OnCell()))))
    #println("interpolator used")
    x = X * (n-1)  .+ 1;
    y = Y * (n-1)  .+ 1;
    return cdf(x,y);
end

function (obj::copula)(X :: Sampleable{Univariate},Y :: Sampleable{Univariate})             #Constructing a joint distribution using the copula
    return joint(X,Y,obj)
end

function sample(C :: AbstractCopula, N = 1; plot = false)

    if plot return samplePlot(C,N);end
    if (C.func == Gau) return CholeskyGaussian(N, C.param) ;end     # Use Cholesky decompostition of the Cov matrix for gaussian copula sampling

    x = rand(N);    y = rand(N);
    ux = x;         uy = zeros(N);
    m = n;          useInterp = false;      # Set true to use interpolator

    if C.func   == perf return hcat(ux,ux); end
    if C.func   == opp return hcat(ux,1 .- ux); end
    if C.func   == indep return hcat(ux,y); end

    # if (!ismissing(C.func) && n < 10^5) m = 10^5; end  # If function is provided more acurate sampling, otherwise must use interpolator
    e = sqrt(eps());      ygrid = range(0,1,length = m);

    #if (C.func == Gau) useInterp = true ;end           # Should you not want to use Cholesky and interpolator has better performance for gaussian

    conditional = (C(x .+ e/2,ygrid, useInterp) - C(x .- e/2,ygrid, useInterp))./e

    for i=1:N
        #conditional = (C(x[i] + e/2,ygrid, true) - C(x[i] - e/2,ygrid, true))./e
        #conditional[1] = 0; conditional[end]=1;       # Small errors sometimes happen
        #if (conditional[end] < y[i]) uy[i] = 1;      # Should linearly interpolate between cond[end] and 1
        #elseif (conditional[1] > y[i]) uy[i] = 0;      # Should linearly interpolate between cond[end] and 1
        #else uy[i] = findlast(conditional .<= y[i]) / m; end        # Finds Nothing??
        #yitp_linear = interpolate(conditional, BSpline(Linear()))
        #uy[i] = yitp_linear(y[i] * m);

        #uy[i] = findlast(conditional .<= y[i])/m; end

        uy[i] = linInterp(conditional[i,:], y[i]);
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

function conditional(C :: AbstractCopula, xVal :: Real; plot = false)   # May also work for joint? xVal will be take from invCdf(M1, u1)

    e = 1/n;        ygrid = range(0,1,length = n);
    if !ismissing(C.func) e = (sqrt(eps())); end

    conditional = (C(xVal + e/2,ygrid) - C(xVal - e/2,ygrid))./e;
    if (plot) plot2(C.cdf, conditional, xVal) end
    return conditional
end

function calcCdf(density :: Array{<:Real,2})
    cdf = zeros(size(density))
    for i = 2:n
        for j = 2:n
            cdf[i,j] =  (density[i,j]/(n^2) + cdf[i, j-1] + cdf[i-1, j] - cdf[i-1, j-1]);
        end
    end
    return cdf
end

function calcDensity(C :: AbstractCopula)

    if ismissing(C.func) return calcDensity(C.cdf); end
    #density = zeros(n,n);
    e = sqrt(sqrt(eps()));
    #x = range(2*e ,1 - 3e    ,length = n);
    #y = range(2*e ,1 - 3e    ,length = n);

    x = range(0 ,1 - e    ,length = n);
    y = range(0 ,1 - e    ,length = n);


    der1 = ( C(x .+ e, y)       - C(x, y)     )/e      # df/dx
    der2 = ( C(x .+ e, y .+ e ) - C(x, y .+ e) )/e     # d(df/dx)/dy
    density = (der2-der1)/e

    return density
end

function calcDensity(cdf :: Array{<:Real,2}, xRange = range(0,1,length = n), yRange = range(0,1,length = n))

    xStep = Float64(xRange.step); yStep = Float64(yRange.step);

    density = zeros(size(cdf));
    for i = 1:n-1
        for j = 1:n-1
            density[i,j] = (cdf[i,j] + cdf[i+1,j+1] - cdf[i+1,j] - cdf[i,j+1])/(xStep*yStep);

            #   -> Possible alternative:
            #der1 = ( cdf[i + 1, j]      - cdf[i, j]     )/(xStep^2)             # df/dx
            #der2 = ( cdf[i + 1, j + 1 ] - cdf[i, j + 1] )/(yStep^2)             # d(df/dx)/dy
            #density[i,j] = (der2 - der1)
        end
    end

    # Problem: How do you get the end values?
    #   Two solutions:

    #   -> If copula, it's symetric (exact answer exept at two endpoints)
    #density[end,:] = reverse(density[:,1]);
    #density[:,end] = reverse(density[1,:]);

    #   -> If joint, it may not be symetric
    density[end,:] = density[end-1,:];
    density[:,end] = density[:,end-1];
    return density
end

function M()
    x = y = range(0,1,length = n);
    return copula(perf(x,y), density = Matrix{Float64}(I, n, n), func = perf);
end

function W()
    x = y = range(0,1,length = n);

    den = zeros(n,n)
    for i = 1:n
        for j = 1:n
            if (i-1) == n-(j-1); den[i,j] = 1; end
        end
    end
    return copula(opp(x,y), density = den , func = opp);
end

function π()
    x = y = range(0,1,length = n);
    return copula(indep(x,y), func = indep);
end
Pi() = π()

function Frank(s = 1)                                     #   Frank copula
    x = y = range(0,1,length = n);                          #   s>0; 0 for perfect, 1 for indep, inf for oposite
    return copula(F(x,y,s), func = F, param = s);
end

function Clayton(t = 0)                                  #   Clayton copula
    x = y = range(0,1,length = n);                          #   t>-1; -1 for opposite, 0 for indep and inf for perfect
    return copula(Cla(x,y,t), func = Cla, param = t);
end

function Gaussian(corr = 0)
    x = y = range(0,1,length = n);
    cdf = Gau(x,y,corr);
    #cdf[:,end] = cdf[end,:] = x;    # This removes the NaNs from Liams mvNormCdf. We know the marginals must be uniform
    return copula(cdf, func = Gau, param = corr);
end

indep(X, Y) = [x*y for x in X, y in Y];
perf(X, Y)  = [min(x,y) for x in X, y in Y];
opp(X, Y)   = [max(x+y-1,0) for x in X, y in Y];
F(X,Y,s = 1)      = [log(1+(s^x-1)*(s^y-1)/(s-1))/log(s) for x in X, y in Y]
Cla(X,Y, t = 0) = [max((x^(-t)+y^(-t)-1)^(-1/t),0) for x in X, y in Y]
Gau(X,Y, corr = 0)  = [bivariate_cdf(quantile.(Normal(),x),quantile.(Normal(),y), corr) for x in X, y in Y];

Arch1(X,Y,theta = 0) = [1 - ((1-x)^theta + (1-y)^theta - (1-x)^theta * (1-y)^theta )^theta for x in X, y in Y];

mvNormCdf(X,Y,cor) = [bivariate_cdf(x,y,cor) for x in X, y in Y];

wait_for_key() = (print(stdout, "press enter to continue"); read(stdin, 1); nothing);

function samplePlot(C :: AbstractCopula, N = 1)

    x = rand(N);    y = rand(N);
    ux = x;         uy = zeros(N);
    m = n;
    e = 1/m;        ygrid = range(0,1,length = m);

    if !ismissing(C.func) e = sqrt(eps()); end

    for i=1:N
        conditional = (C(x[i] + e/2,ygrid) - C(x[i] - e/2,ygrid))./e
        uy[i] = findlast(conditional .<= y[i]) / m;
        println("-------------------------------------------")
        println("N: $i.")
        println("Independent ~U(0,1)'s:")
        println(hcat(x[i],y[i]))
        println();

        plot2(C.cdf,conditional,ux[i], uy[i], y[i])

        println("Samples from copula:");
        println(hcat(ux[i],uy[i]))
        println("-------------------------------------------")

        wait_for_key()
    end
    return hcat(ux,uy);
end

plotDensity(x :: AbstractCopula) = plot1(x.density);
plotCdf(x::AbstractCopula) = plot1(x.cdf);

function plot1(A :: Array{<:Real,2})

#   https://github.com/gizmaa/Julia_Examples/blob/master/pyplot_surfaceplot.jl
    m = size(A)[1];
    if m < pn; ppn = m; else ppn = pn; end

    x = y = range(0,1,length=ppn)
    xgrid = repeat(x',ppn,1)
    ygrid = repeat(y,1,ppn)

    nm = round(m/ppn);

    z = A[1:Int(nm):end,1:Int(nm):end]

    fig = figure("SurfacPlots",figsize=(10,10))
    ax = fig.add_subplot(2,1,1,projection="3d")
    plot_surface(xgrid, ygrid, z, rstride=2,edgecolors="k", cstride=2, alpha=0.8, linewidth=0.25,cmap=ColorMap("coolwarm"))
    xlabel("X")
    ylabel("Y")
    PyPlot.title("Surface Plot")

    subplot(212)
    ax = fig.add_subplot(2,1,2)
    cp = contour(xgrid, ygrid, z,cmap=ColorMap("coolwarm"),levels = 15)
    ax.clabel(cp, inline=1, fontsize=10)
    xlabel("X")
    ylabel("Y")
    PyPlot.title("Contour Plot")
    tight_layout()

end

function plot2(A :: Array{<:Real,2}, cond:: Array{<:Real,1}, ux = missing, uy = missing, uuy = missing)

    if !isempty(get_fignums()) PyPlot.clf(); end

    m = size(A)[1];
    if m < pn; ppn = m; else ppn = pn; end

    x = y = range(0,1,length=ppn)
    xgrid = repeat(x',ppn,1)
    ygrid = repeat(y,1,ppn)

    nm = m/ppn;

    z = A[1:Int(nm):end, 1:Int(nm):end]
    conditional = cond[1:Int(nm):end];

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

function plot3(A :: Array{<:Real,2}, x = range(0,1,length = pn), y = range(0,1,length = pn))    # Won't work if n!=200

#   https://github.com/gizmaa/Julia_Examples/blob/master/pyplot_surfaceplot.jl

    #m1 = size(A)[1];
    #m2 = size(A)[2];
    #if m1 < pn; ppn1 = m1; else ppn1 = pn; end      # ppn1 is the number of points to plot in X
    #if m2 < pn; ppn2 = m2; else ppn2 = pn; end      # ppn2 is the number of points to plot in Y

    #nm = round(m/ppn);                              # The descritization of the domain

    #x = range(rangeX[1],rangeX[2],length=pn)
    #y = range(rangeY[1],rangeY[1],length=pn)
    xgrid = repeat(x',pn,1)
    ygrid = repeat(y,1,pn)
    z = A'

    fig = figure("SurfacePlots",figsize=(10,10))
    ax = fig.add_subplot(2,1,1,projection="3d")
    plot_surface(xgrid, ygrid, z, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("coolwarm"), alpha=0.8, linewidth=0.25)
    #plot_surface(xgrid, ygrid, z, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("Gray"), alpha=0.8, linewidth=0.25)
    xlabel("X")
    ylabel("Y")
    PyPlot.title("Surface Plot")

    #subplot(212)
    ax = fig.add_subplot(2,1,2)
    cp = contour(xgrid, ygrid, z, cmap=ColorMap("coolwarm"))
    ax.clabel(cp, inline=1, fontsize=10)
    xlabel("X")
    ylabel("Y")
    PyPlot.title("Contour Plot")
    tight_layout()

end

function plot4(J :: AbstractJoint, CDF = true)    # Won't work if n!=200

#   https://github.com/gizmaa/Julia_Examples/blob/master/pyplot_surfaceplot.jl

    x = J.xRange;       y = J.yRange;
    M1 = J.marginal1;   M2 = J.marginal2;
    xgrid = repeat(x',pn,1)
    ygrid = repeat(y,1,pn)
    if CDF; z = J.cdf; title = "cdf" else z = J.density; title = "pdf";end


    fig1 = figure("SurfacePlots",figsize=(10,10))
    #ax = fig.add_subplot(2,1,1,projection="3d")
    plot_surface(xgrid, ygrid, z, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("coolwarm"), alpha=0.8, linewidth=0.25)
    xlabel("X")
    ylabel("Y")

    PyPlot.title("$title surface plot")

    fig = plt.figure("ContourPlots", figsize=(10, 10))
    grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
    main_ax = fig.add_subplot(get(grid, (slice(1,4),slice(0,3))))
    yMarg  = fig.add_subplot(get(grid, (slice(1,4),3)), xticklabels=[], sharey=main_ax)
    xMarg  = fig.add_subplot(get(grid, (0,slice(0,3))), yticklabels=[], sharex=main_ax)

    cp = main_ax.contour(xgrid, ygrid, z, cmap=ColorMap("coolwarm"), levels =15)
    main_ax.clabel(cp, inline=1, fontsize=10)
    #xlabel("X")
    #ylabel("Y")
    if (CDF) xMarg.plot(x,cdf.(M1,x));else xMarg.plot(x,pdf.(M1,x));end
    if (CDF) yMarg.plot(cdf.(M2,y),y);else yMarg.plot(pdf.(M2,y),y);end

    PyPlot.title("joint $title and marginals")
    tight_layout()

end

slice(x,y) = pycall(pybuiltin("slice"), PyObject, x, y)         # For performing python array slicing for

function scatter(a :: Array{Float64,2}; title = "samples")

    x = a[:,1];
    y = a[:,2];

    fig = plt.figure(title, figsize=(10, 10))
    grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
    main_ax = fig.add_subplot(get(grid, (slice(1,4),slice(0,3))))
    y_hist  = fig.add_subplot(get(grid, (slice(1,4),3)), xticklabels=[], sharey=main_ax)
    x_hist  = fig.add_subplot(get(grid, (0,slice(0,3))), yticklabels=[], sharex=main_ax)

    # scatter points on the main axes
    main_ax.plot(x, y, "o", markersize=3, alpha=0.2)

    # histogram on the attached axes
    x_hist.hist(x, 40, histtype="stepfilled", orientation="vertical", color="gray")
    #x_hist.invert_yaxis()

    y_hist.hist(y, 40, histtype="stepfilled", orientation="horizontal", color="gray")
    #y_hist.invert_xaxis()

end



mutable struct joint <: AbstractJoint

    marginal1 :: Sampleable{Univariate}
    marginal2 :: Sampleable{Univariate}
    copula  :: AbstractCopula

    xRange
    yRange

    cdf ::  Union{Array{Float64,2},Missing}     # Or should the density/ cdf be calculated on the fly from the copula?
    density :: Union{Array{Float64,2},Missing}

    function joint(marginal1, marginal2, copula)

        #if (ismissing(cdf) && ismissing(density)) throw(ArgumentError("Cdf and/or density must be provided")) end
        #if (ismissing(cdf)); cdf = calcCdf(density); end
        #if (ismissing(density)); density = calcDensity(cdf);end

        xRange = range(invCdf(marginal1,bOt), invCdf(marginal1,tOp), length = n)
        yRange = range(invCdf(marginal2,bOt), invCdf(marginal2,tOp), length = n)

        j = new(marginal1, marginal2, copula, xRange, yRange)
        cdf = j(xRange,yRange)
        j = new(marginal1, marginal2, copula, xRange, yRange, cdf)
        density = calcDensity(j)

        return new(marginal1, marginal2, copula, xRange, yRange, cdf, density)

    end
end

function (obj :: joint)(X,Y)    # Evaluating the joint cdf
                                # BiLinear interpolator if cdf is provided
    return obj.copula(cdf.(obj.marginal1,X),cdf.(obj.marginal2,Y))  # This is the case where we have the functions
end

function sample(J :: AbstractJoint, N =1)

    copulaSamples = sample(J.copula,N)

    x = invCdf(J.marginal1,copulaSamples[:,1]);
    y = invCdf(J.marginal2,copulaSamples[:,2]);
    return hcat(x,y)

end

function invCdf(X :: Sampleable{Univariate}, u)
    return quantile.(X,u)                           # We will also need the interpolartor here
end

plotCdf(j::AbstractJoint) = plot3(j.cdf, j.xRange,j.yRange)
plotDensity(j::AbstractJoint) = plot3(j.density, j.xRange,j.yRange)


function calcDensity(J :: AbstractJoint)        # Legacy script. All it does not is return calcDensity(j.cdf,xRange,yRange);

    e = sqrt(sqrt(eps()));

    x = J.xRange .- e;
    y = J.yRange .- e;

    der1 = ( J(x .+ e, y)       - J(x, y)     )/e     # df/dx
    der2 = ( J(x .+ e, y .+ e ) - J(x, y .+ e) )/e     # d(df/dx)/dy
    density = (der2-der1)/e

    return density

end

M(M1::Sampleable{Univariate}, M2::Sampleable{Univariate}) = M()(M1,M2)
W(M1::Sampleable{Univariate}, M2::Sampleable{Univariate}) = W()(M1,M2)
π(M1::Sampleable{Univariate}, M2::Sampleable{Univariate}) = π()(M1,M2)
Frank(M1::Sampleable{Univariate}, M2::Sampleable{Univariate}, s = 1) = Frank(s)(M1,M2)
Clayton(M1::Sampleable{Univariate}, M2::Sampleable{Univariate}, t = 0) = Clayton(t)(M1,M2)
Gaussian(M1::Sampleable{Univariate}, M2::Sampleable{Univariate}, corr = 0) = Gaussian(corr)(M1,M2)

function Base.show(io::IO, z::copula)

    statement1 = "Arbitrary"
    statement2 = ""

    if (!ismissing(z.func))
        func = z.func
        if (func == indep); func ="π";end
        if (func == perf); func = "M";end
        if (func == opp); func = "W";end
        statement1 = "$(func)"
    end

    if (!ismissing(z.param))
        func = z.func
        parName = "par"
        if (func == Gau) parName = "r";end
        if (func == F) parName = "s";end
        if (func == Cla) parName = "t";end
        statement2 = "$parName=$(z.param)"
    end

    print(io, "Copula ~ $statement1($statement2)");
end

function Base.show(io::IO, z::joint)

    statement1 = "Arbitrary"
    statement2 = ""

    if (!ismissing(z.copula.func))
        func = z.copula.func
        if (func == indep); func ="π";end
        if (func == perf); func = "M";end
        if (func == opp); func = "W";end
        statement1 = "$(func)"
    end

    if (!ismissing(z.copula.param))
        func = z.copula.func
        parName = "par"
        if (func == Gau) parName = "r";end
        if (func == F) parName = "s";end
        if (func == Cla) parName = "t";end
        statement2 = "$parName=$(z.copula.param), "
    end

    print(io, "Joint ~ $statement1( $statement2$(z.marginal1), $(z.marginal2) )");
end

function linInterp(A, x)

    ind0 = findlast(A .<= x);
    x0 = A[ind0];   x1 = A[ind0+1]
    y0 = ind0/n;     y1 = (ind0+1)/n

    return y0 + (x-x0) * (y1-y0)/(x1-x0)
end
#=
function calcDensity(J :: AbstractJoint)        # Legacy script. All it does not is return calcDensity(j.cdf,xRange,yRange);

    #return calcDensity(J.cdf,J.xRange,J.yRange);

    if (typeof(C) <: AbstractJoint)
        func = C.copula.func;
        xRange = C.xRange; yRange = C.yRange;
    else
        func = C.func;
        xRange = yRange = range(0,1,length = n);
    end

    if ismissing(func) return calcDensity(C.cdf, xRange, yRange); end
    density = zeros(n,n);
    e = sqrt(sqrt(eps()));
    x = range(2*e + Float64(xRange.ref), Float64(xRange.offset) - 3e    ,length = xRange.len);
    y = range(2*e + Float64(yRange.ref), Float64(yRange.offset) - 3e    ,length = yRange.len);
    for i = 1:n
        for j = 1:n
            der1 = ( C(x[i] + e, y[j])      - C(x[i], y[j])     )/e     # df/dx
            der2 = ( C(x[i] + e, y[j] + e ) - C(x[i], y[j] + e) )/e     # d(df/dx)/dy
            density[i,j] = (der2-der1)/e
        end
    end
    return density

end
=#
#=
function bilinearInterp(A :: Array{<:Real,2}, x:: Array{<:Real,1})

    if (length(x) != ndims(A)) throw(ArgumentError("Dimension of function and point must both be 2")); end
    if (length(x) != 2)        throw(ArgumentError("Dimension point must both be 2")); end
    if (ndims(A) != 2)        throw(ArgumentError("Dimension function must both be 2")); end

end

=#

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
