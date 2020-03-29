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
#   Port of R code pba.r by Scott Ferson and Jason O'Rawe, Applied Biomathematics
#   Origional code available at: https://github.com/ScottFerson/pba.r
######


#####
#   Known Bugs:
#
#   ->  After using a pbox in a computation, the shape is no loger saved. For example normal(0,1) + 1, should still be a normal
#   ->  env() between two scalars should return interval. Works but returned variance is [0,0]
#   ->  Giving pboxes to the interval constructor should return envelope
#   ->  Potential error with reciprocate function. No mean and variance transformations yet
#   ->  Cut(pbox([0,1]),0.5) returns [0.5026256564141035, 0.4973743435858965]. Right side larger!
#   ->  Multiplication of pboxes with intervals not working as expected.
#
#   Severe:
#       ->  makepbox() won't work when both a pbox and an interval type are introduced as arguments
###

#####
#   Fixed Bugs wrt. pba.r:
#
#   ->  When passing a pbox to the pbox constructor, mean and variance not saved and recalculated from bounds
#
###

abstract type AbstractPbox <: Real end

mutable struct pbox <: AbstractPbox
    id :: String                        # Id of the pbox ie: "PB 1"
    u :: Union{Array{<:Real},Real}      # vector of the left bounds
    d :: Union{Array{<:Real},Real}      # vector of the right bounds
    n :: Int64                          # descretization size
    shape :: String                     # shape of internal distribtion
    ml :: Real                          # lower bound of mean
    mh :: Real                          # upper bound of mean
    vl :: Real                          # lower bound on variance
    vh :: Real                          # upper bound of variance
    name :: String                      # name
    dids :: String
    bob :: Int64                        # for dependancy tracking


    function pbox(u::Union{Missing, Array{<:Real}, Real} = missing, d=u; shape = "", name = "", ml= missing, mh = missing,
        vl = missing, vh = missing, interpolation = "linear",
        bob = missing, perfect = missing, opposite = missing, depends = missing, dids=missing)

        steps = ProbabilityBoundsAnalysis.steps;                  # Reading global varaibles costly, creating local for multiple use

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
        else

            if (!(typeof(u)<:Union{Array{<:Real},Real}) || !(typeof(d)<:Union{Array{<:Real},Real}))
                throw(ArgumentError("Bounds of a p-box must be numeric"))
            end

            if (typeof(u)<:AbstractInterval) u = u.lo; end
            if (typeof(d)<:AbstractInterval) d = d.hi; end

            # Should we always be interpolating? What if we convolve 2 pboxes of different steps and ProbabilityBoundsAnalysis.steps is larger
            # The resulting pbox should have the smaller of the 3
            if (length(u) != steps) u = Interpolate(u,interpolation,true);end
            if (length(d) != steps) d = Interpolate(d,interpolation,false);end

            unique = uniquePbox();
            id = "PB $unique";

            p = new(id,u,d,steps,"",-∞,∞,0,∞,"","",unique);
        end

        p = computeMoments(p);

        if (!ismissing(shape))      p.shape = shape; end
        if (!ismissing(ml))         p.ml = max(ml,p.ml);end
        if (!ismissing(mh))         p.mh = min(mh,p.mh);end
        if (!ismissing(vl))         p.vl = max(vl,p.vl);end
        if (!ismissing(vh))         p.vh = min(vh,p.vh);end
        if (!ismissing(dids))       p.dids = "$(p.id) $dids";end
        if (!ismissing(opposite))   p.bob = -opposite.bob;end
        if (!ismissing(perfect))    p.bob = perfect.bob;end
        if (!(ismissing(depends)))  p.dids = "$(p.id) $(depends.dids)";end
        if (!ismissing(name))       p.name = name;end
        #if (ProbabilityBoundsAnalysis.plottingEvery)      plotpbox(p); end

        p = checkMoments(p);

        return p;
    end

end

# Outer constructor

#=
pbox(;u::Union{Missing, Array{<:Real}, Real} = missing, d=u, shape = "", name = "", ml= missing, mh = missing,
    vl = missing, vh = missing, interpolation = "linear",
    bob = missing, perfect = missing, opposite = missing, depends = missing, dids=missing) =
    pbox(u,d,shape = shape, name = name, ml = ml, mh = mh,
    vl = vl, vh = vh, interpolation = interpolation,
    bob = bob, perfect = perfect, opposite = opposite, depends = depends, dids = dids);

=#
uniquePbox() = ( ProbabilityBoundsAnalysis.setPboxes(ProbabilityBoundsAnalysis.pboxes+1); return ProbabilityBoundsAnalysis.pboxes );

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

#pbox(x :: pbox) = pbox(x);


#sideVariance(w :: Array{<:Real}, mu=missing) = if (ismissing(mu)) mu = mean(w);end return (max(0,mean((w .- mu).^2)));end


####################################
# Pbox moments                     #
####################################

mean(x::pbox) = Interval(x.ml,x.mh);
var(x::pbox) = Interval(x.vl,x.vh);
std(x::pbox) = √(var(x));

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
      if (1<ProbabilityBoundsAnalysis.verbose) println(ProbabilityBoundsAnalysis.meanDisagreementMessage);end
  end

  a = var(x);
  b = dwVariance(x);
  x.vl = max(left(a),left(b))
  x.vh = min(right(a),right(b))

  if (x.vh < x.vl)
    # use the observed mean
    x.vl = left(b)
    x.vh = right(b)
    if (1<ProbabilityBoundsAnalysis.verbose) println(ProbabilityBoundsAnalysis.varDisagreementMessage);end
end

    return x;

end


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

function Interpolate(x::Union{Array{<:Real},Real}, interpolation::String = "linear", left::Bool = true)

    if (interpolation == "linear") return linearInterpolation(x);
    elseif (interpolation == "outer") return outerInterpolation(x,left);
    elseif (interpolation == "step") return stepInterpolation(x);
    elseif (interpolation == "spline") return splineInterpolation(x);
    else throw(ArgumentError("Interpolation option not valid. The available schemes are linear, outer, step and spline"));
    end

end

function linearInterpolation( V::Union{Array{<:Real},Real} )

    steps = ProbabilityBoundsAnalysis.steps;

    m = length(V) - 1;
    if (m == 0) return ([V for i=1:steps]);end
    if (ProbabilityBoundsAnalysis.steps == 1) return [min(V), max(V)];end

    d = 1/m;
    n = Int(round(d * steps * 20));
    if (n==0) c = V;
    else
        c = [];
        for i = 1:m
            v = V[i];
            w = V[i+1];
            c = [c;LinRange(v,w,n)];
        end
    end
    u = zeros(1,steps);
    for k = 1:steps; u[k] = c[1+Int(round((length(c)-1)*(k-1)/(steps-1)))]; end
    return u;
end

###
#   Still needed
###

#   outerInterpolation
#   stepInterpolation
#   splineInterpolation
