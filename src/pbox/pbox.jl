######
# This file is part of the ProbabilityBoundsAnalysis.jl package
#
#   Definition of pbox class with constructor methdos
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
#
#   About 50% of this file is a port of R code pba.r by Scott Ferson and Jason O'Rawe, Applied Biomathematics
#   Origional code available at: https://github.com/ScottFerson/pba.r
######

######
#   Known Bugs:
#
#   ->  After using a pbox in a computation, the shape is no longer saved. For example normal(0,1) + 1, should still be a normal
#   ->  env() between two scalars should return interval. Works but returned variance is [0,0]
#   ->  Giving pboxes to the interval constructor should return envelope
#   ->  Potential error with reciprocate function. No mean and variance transformations yet
#   ->  Cut(pbox([0,1]),0.5) returns [0.5026256564141035, 0.4973743435858965]. Right side larger!
#   ->  Multiplication of pboxes with intervals not working as expected.
#   ->  id and dids don't work. When a pbox is made, 4 is added to the id counter
#
#   ->  Check bounded. What is the reciprocal of bounded? What is it's complement?
#       -> Does my arithmetic with bounded make sense?
#
#   Severe:
#       ->  makepbox() won't work when both a pbox and an interval type are introduced as arguments
#####

#####
#   Fixed Bugs wrt. pba.r:
#
#   ->  When passing a pbox to the pbox constructor, mean and variance not saved and recalculated from bounds
#
#####


"""
    mutable struct pbox <: AbstractPbox

Basic type used in Probability Bounds Analysis. A set of probability distributions defined by an upper (u) and lower (d/down) cdf, interval bounds on mean and variance and distibution shape (eg. Normal).

If partial information is given (eg, no bounds on moments), the other properties are determined from provided information.

# Constructors
* `pbox()                     => Vacous case`
* `pbox(u :: Real)            => Scalar case`
* `pbox(u :: Real, d :: Real) => Interval case`
* `pbox(u :: Array{<:Real}, d :: Array{<:Real}, shape :: String, ml :: Real, mh :: Real, vl :: Real, vh :: Real)`

where ml/mh are the lower/upper mean, and vl/vh are the lower/upper variance.

See also: [`mean`](@ref), [`var`](@ref), [`uniform`](@ref), [`normal`](@ref), [`conv`](@ref), [`copula`](@ref)
"""
mutable struct pbox <: AbstractPbox
    u :: Union{Array{<:Real},Real}      # vector of the left bounds
    d :: Union{Array{<:Real},Real}      # vector of the right bounds
    n :: Int64                          # descretization size
    shape :: String                     # shape of internal distribtion
    ml :: Real                          # lower bound of mean
    mh :: Real                          # upper bound of mean
    vl :: Real                          # lower bound on variance
    vh :: Real                          # upper bound of variance
    name :: String                      # name
    bounded :: Array{Bool,1}            # Is it bounded?

    function pbox(u::Union{Missing, Array{<:Real}, Real} = missing, d=u; shape = "", name = "", ml= missing, mh = missing,
        vl = missing, vh = missing, bounded = [false, false])

        steps = parametersPBA.steps;                  # Reading global varaibles costly, creating local for multiple use

        if (ismissing(u) && ismissing(d))
            u = -∞;
            d = ∞;
        end

        if (typeof(u) == pbox)
            p = u;
            ml = p.ml;
            mh = p.mh;
            vl = p.vl;
            vh = p.vh;
            shape = p.shape;
            bounded = p.bounded;
        else

            if (!(typeof(u)<:Union{Array{<:Real},Real}) || !(typeof(d)<:Union{Array{<:Real},Real}))
                throw(ArgumentError("Bounds of a p-box must be numeric"))
            end

            if (typeof(u)<:AbstractInterval) u = u.lo; end
            if (typeof(d)<:AbstractInterval) d = d.hi; end

            iis = ProbabilityBoundsAnalysis.iii()
            jjs = ProbabilityBoundsAnalysis.jjj()

            if bounded[1]; iis = ProbabilityBoundsAnalysis.ii(); end
            if bounded[2]; jjs = ProbabilityBoundsAnalysis.jj(); end

            # Should we always be interpolating? What if we convolve 2 pboxes of different steps and ProbabilityBoundsAnalysis.steps is larger
            # The resulting pbox should have the smaller of the 3
            if (length(u) != steps) u = linearInterpolation(u, iis); end
            if (length(d) != steps) d = linearInterpolation(d, jjs); end

            p = new(u,d,steps,"",-∞,∞,0,∞,"", bounded);
        end

        p = computeMoments(p);

        if (!ismissing(shape))      p.shape = shape; end
        if (!ismissing(ml))         p.ml = max(ml,p.ml);end
        if (!ismissing(mh))         p.mh = min(mh,p.mh);end
        if (!ismissing(vl))         p.vl = max(vl,p.vl);end
        if (!ismissing(vh))         p.vh = min(vh,p.vh);end
        if (!ismissing(name))       p.name = name;end

        p = checkMoments(p);

        return p;
    end

end

###
#   If type is called, return cdf
###
function (obj::pbox)(x)
    return cdf(obj,x)
end

###
# Returns bounds on cdf. Upper bound may not be best possible
###
"""
	cdf(s :: pbox, x :: Real)

Returns cdf value (interval) at x

See also: [`cut`](@ref), [`mass`](@ref), [`rand`](@ref)
"""
function cdf(s :: pbox, x::Real)
    d = s.d; u = s.u; n = s.n;
    bounded = s.bounded;

    if x < u[1];
        return interval(0,1/n) * (1-bounded[1]);
    end;
    if x > d[end];
        if bounded[2]; return 1; end
        return interval((n-1)/n, 1);
    end;

    indUb = 1 - sum(x .< u)/n;
    indLb = 1 - sum(x .<= d)/n;

    pub = 1; plb = 0;

    #if x < u[end];  pub = indUb; end
    #if x > d[1];    plb = indLb; end

    #return interval(plb, pub)

    return interval(indLb, indUb)
end


"""
	cdf(s :: pbox, x :: Interval)

Returns interval bounds on cdf value in interval x

See also: [`cut`](@ref), [`mass`](@ref), [`rand`](@ref)
"""
function cdf(s :: pbox, x::Interval)
    lb = left(cdf(s,left(x)));
    ub = right(cdf(s,right(x)));
    return interval(lb, ub)
end


"""
	mass(s :: pbox, lo :: Real, hi :: Real)

Returns bounds on probability mass in interval [lo, hi]

See also: [`cut`](@ref), [`cdf`](@ref), [`rand`](@ref)
"""
function mass(s :: pbox, lo :: Real, hi :: Real)

    cdfHi = cdf(s, hi)
    cdfLo = cdf(s, lo)

    ub = min(1, right(cdfHi) - left(cdfLo))
    lb = max(0, left(cdfHi) - right(cdfLo))

    return interval(lb, ub)

end

"""
	mass(s :: pbox, x:: Interval)

Returns bounds on probability mass in interval x

See also: [`cut`](@ref), [`cdf`](@ref), [`rand`](@ref)
"""
mass(s :: pbox, x:: Interval) = mass(s, x.lo, x.hi)


"""
    makepbox(x...)

Returns an array of pboxes from an array of inputs (eg an array of intervals or reals).

# Examples
```jldoctest
julia> s = makepbox(interval(0,1))
Pbox: 	  ~  ( range=[0.0, 1.0], mean=[0.0, 1.0], var=[0.0, 0.25])

julia> array = [interval(0, 1), interval(0, 2), 3];

julia> s = makepbox.(array)
3-element Vector{pbox}:
 Pbox: 	  ~  ( range=[0.0, 1.0], mean=[0.0, 1.0], var=[0.0, 0.25])
 Pbox: 	  ~  ( range=[0.0, 2.0], mean=[0.0, 2.0], var=[0.0, 1.0])
 Pbox: 	  ~  ( range=3.0, mean=3.0, var=0.0)
```
"""
function makepbox(x...)


    # This line is causing alot of problems with the interval package. Interval package converts all elements of the array to the same type Interval{<:Float64} etc
    elts :: Array{Any} = [x...];
    # Not sure about this line, if you input a Tuple you can still use the function as expected
    # if (typeof(x[1])<:Tuple) elts = x[1]; end
    if (length(elts)==1) return pbox(x...);end
    # elts = convert(Array{Any},elts);
    for i=1:length(elts)
        if (!ispbox(elts[i]))
            elts[i] = pbox(elts[i]);
        end
    end
    return elts;
end


"""
    pbox( x :: Array{Interval{T}, 1} ) where T <: Real

Constructs a pbox from an array of intervals with equal mass. Left and right bounds are sorted to construct cdf bounds.
"""
function pbox( x :: Array{Interval{T}, 1}, bounded :: Vector{Bool} = [true, true]) where T <: Real

    us = left.(x);  ds = right.(x)
    us = sort(us);  ds = sort(ds)

    numel = length(us)

    n = parametersPBA.steps;

    uNew = zeros(n);    dNew = zeros(n);

    i = range(0,stop = 1, length = numel +1);
    j = i[1:end-1];   i = i[2:end];

    j = reverse(j)

    iis = ProbabilityBoundsAnalysis.ii();
    jjs = ProbabilityBoundsAnalysis.jj();

    #jjs = reverse(jjs)

    for k = 1:n

        iThis = findfirst(iis[k] .<= i)
        jThis = findfirst(jjs[k] .>= j)

        uNew[k] = us[iThis]
        dNew[k] = ds[jThis]

    end

    dNew = reverse(dNew)

    return pbox(uNew, dNew, bounded = bounded)

end


##
#   Making a pbox from a matrix of 2nd order samples [Nouter x Ninner]
##
"""
    pbox( x :: Array{T, 2} ) where T <: Real

Constructs a pbox from a matrix of 2nd order samples [Nouter x Ninner]
"""
function pbox( x :: Array{T, 2}, bounded :: Vector{Bool} = [true, true]) where T <: Real

    x = sort(x, dims = 2);

    Nouter, Ninner = size(x);
    us = minimum(x, dims = 1);
    ds = maximum(x, dims = 1);

    us = sort(us', dims=1)
    ds = sort(ds', dims=1)

    n = parametersPBA.steps;

    uNew = zeros(n);    dNew = zeros(n);

    i = range(0,stop = 1, length = Ninner + 1);
    j = i[1:end-1];   i = i[2:end];

    j = reverse(j)

    iis = ProbabilityBoundsAnalysis.ii();
    jjs = ProbabilityBoundsAnalysis.jj();

    #jjs = reverse(jjs)

    for k = 1:n

        iThis = findfirst(iis[k] .<= i)
        jThis = findfirst(jjs[k] .>= j)

        uNew[k] = us[iThis]
        dNew[k] = ds[jThis]

    end

    dNew = reverse(dNew)

    return pbox(uNew, dNew, bounded = bounded)

end

#pbox(x :: pbox) = pbox(x);


#sideVariance(w :: Array{<:Real}, mu=missing) = if (ismissing(mu)) mu = mean(w);end return (max(0,mean((w .- mu).^2)));end


####################################
# Pbox moments                     #
####################################

"""
    mean( x :: pbox)

Get interval mean of a pbox
See also: [`var`](@ref), [`std`](@ref)
"""
mean(x::pbox) = Interval(x.ml,x.mh);

"""
    var( x :: pbox)

Get interval variance of a pbox
See also: [`std`](@ref), [`mean`](@ref)
"""
var(x::pbox) = Interval(x.vl,x.vh);

"""
    std( x :: pbox)

Get interval std of a pbox
See also: [`var`](@ref), [`mean`](@ref)
"""
std(x::pbox) = √(var(x));

var(x :: Interval) = var(pbox(x))

dwmean( x :: pbox) = Interval(mean(x.u),mean(x.d));

function dwVariance(x :: pbox)
    if (any(x.u .== -Inf) || any(x.u .== Inf)) return interval(0,∞);end
    if (all(x.u .== x.u[1]) && all(x.d .== x.d[1])) return interval(0,(x.d[1]-x.u[1])^2/4); end

    vr = sideVariance(x.u);
    w = deepcopy(x.u);
    n = length(x.u);
    for i = n:-1:1
        w[i] = x.d[i];
        v = sideVariance(w,mean(w));
        if (isnan(vr) || isnan(v)) vr = Inf elseif (vr<v) vr = v; end
    end

    if (x.u[n] <= x.d[1])
        vl = 0.0;
    else
        w = deepcopy(x.d);
        vl = sideVariance(w,mean(w));
        for i = n:-1:1
             w[i] = x.u[i];
             here = w[i];
             if (1<i)
                 for j = (i-1):-1:1
                     if (w[i] < w[j])
                         w[j] = here
                     end
                 end
             end
             v = sideVariance(w,mean(w))
             if (isnan(vl) || isnan(v)) vl = 0 elseif (v<vl) vl = v;end
        end
    end
    return interval(vl,vr);
end


function sideVariance(w :: Array{<:Real}, mu=missing)
     if (ismissing(mu))
         mu = mean(w);
     end
     return (max(0,mean((w .- mu).^2)));
end



function computeMoments(x :: pbox)

    x.ml = max(x.ml,mean(x.u));
    x.mh = min(x.mh,mean(x.d));
    if (isinterval(x))
        x.vl = max(x.vl,0);
        x.vh = min(x.vh, ((x.u-x.d).^2/4)...); # Should be changed (µ - lo)(µ - hi) if mu is included? What if it is bounded?
        return x;
    end

    if (any(x.u .<= -Inf) || any(x.d .>= Inf)) return x; end    # If range is unbounded don't bother computing variance

    V = 0; JJ = 0;
    j = 1:x.n;
    for J = 0:x.n
        ud = [x.u[j .< J]; x.d[J .<= j]];
        v = sideVariance(ud);
        if ( V < v )
            JJ = J;
            V = v;
        end
    end
    x.vh = V;
    return x;
end

# Check Moments alters the pbox. Also the lower bound of the variance is defined here
# Should this not be done in compute moments?
function checkMoments( x :: pbox)

    a = mean(x);
    b = dwmean(x);
    x.ml = max(left(a),left(b));
    x.mh = min(right(a),right(b));

    if (x.mh < x.ml)
      # use the observed mean
      x.ml = left(b)
      x.mh = right(b)
      if (1<parametersPBA.verbose) println("Disagreement between theoretical and observed mean\n");end
  end

  a = var(x);
  b = dwVariance(x);
  x.vl = max(left(a),left(b))
  x.vh = min(right(a),right(b))

  if (x.vh < x.vl)
    # use the observed mean
    x.vl = left(b)
    x.vh = right(b)
    if (1<parametersPBA.verbose) println("Disagreement between theoretical and observed variance\n");end
end

    return x;

end


######
#   Make stochastic mixture of p-boxes
######

function mixture( x :: Vector{pbox}, w :: Vector{<:Real} = ones(length(x)))

    k = length(x);
    if k != length(w); throw(ArgumentError("Number of weights and mixture elements not equal")); end

    w = w ./ sum(w);

    u, d, n = Float64[], Float64[], Float64[]
    ml, mh, m, vl, vh, v = [], [], [], [], [], []

    for i = 1:k
        u = [u; x[i].u]
        d = [d; x[i].d]
        steps = length(x[i].d)
        n = [n; w[i] * repeat([1/steps], steps)]

        mu = mean(x[i]);
        ml = [ml; left(mu)];
        mh = [mh; right(mu)];
        m  = [m; mu]

        σ2 = var(x[i])
        vl = [vl; left(σ2)]
        vh = [vh; right(σ2)]
        v = [v; σ2]

    end

    n = n / sum(n)
    su = sort(u); su = [su[1]; su]
    pu = [0; cumsum(n[1:length(u)])]

    sd = sort(d); sd = [sd; sd[end]]
    pd = [cumsum(n[1:length(d)]); 1]

    u, d = Float64[], Float64[]
    j = length(pu);

    for p in reverse(ProbabilityBoundsAnalysis.ii())
        while true
            if pu[j] <= p
                break
            end
            j = j - 1
        end
        u = [su[j]; u]
    end

    j = 1

    for p in ProbabilityBoundsAnalysis.jj()
        while true
            if  p <= pd[j]
                break
            end
            j = j + 1
        end
        d = [d; sd[j]]
    end

    mu = interval(sum(w .* ml), sum(w .* mh))

    s2 = 0
    for i = 1:k; s2 = s2 + w[i] * (v[i] + pow(m[i],2)); end
    s2 = s2 - pow(mu, 2)

    return pbox(u, d, ml = left(mu), mh = right(mu), vl = left(s2), vh = right(s2))

end

###
#   Stochastic mixture of intervals
###

function mixture( x :: Vector{Interval{T}}, w :: Vector{<:Real}) where T <: Real

    ws = w ./sum(w);

    lefts = left.(x)
    rights = right.(x)

    ls = sortperm(lefts)
    rs = sortperm(rights)

    lefts = lefts[ls]
    rights = rights[rs]

    wL = ws[ls]
    wR = ws[rs]

    su = lefts; su = [su[1]; su]
    pu = [0; cumsum(wL)]

    sd = rights; sd = [sd; sd[end]]
    pd = [cumsum(wR); 1]

    u, d = Float64[], Float64[]
    j = length(pu);

    for p in reverse(ProbabilityBoundsAnalysis.ii())
        while true
            if pu[j] <= p
                break
            end
            j = j - 1
        end
        u = [su[j]; u]
    end

    j = 1

    for p in ProbabilityBoundsAnalysis.jj()
        while true
            if  p <= pd[j]
                break
            end
            j = j + 1
        end
        d = [d; sd[j]]
    end

    return pbox(u, d, bounded = [true, true])


end

pbox(x :: Vector{Interval{T}}, w :: Vector{<:Real} ) where T <: Real = mixture(x, w)

####################################
# Interpolation schemes            #
####################################

# The pbox() constructor accepts lists of x-values for
# the left and right bounds.  Four different interpolation
# schemes can be used with these lists, including linear,
# spline, step and outward interpolations. The spline
# interpolation scheme requires the interpolation.jl package to be
# loaded.  See also the quantiles() constructor when
# you have quantiles (or percentiles or fractiles), or
# the pointlist() constructor when you have a point
# list description of the p-box.


function linearInterpolation( V::Union{Array{<:Real}, Real}, ps )

    steps = parametersPBA.steps;

    if (length(V) == 1) return ([V for i=1:steps]);end
    if (parametersPBA.steps == 1) return [min(V), max(V)];end

    ipt = interpolate((range(0, 1, length = length(V)),), V, Gridded(Linear()))
    return ipt(ps);
end

###
#   Still needed
###

#   outerInterpolation
#   stepInterpolation
#   splineInterpolation
