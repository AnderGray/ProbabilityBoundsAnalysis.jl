######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Functions to construct pboxes using distributions
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
#
#   About 50% of this file is a port of R code pba.r by Scott Ferson and Jason O'Rawe, Applied Biomathematics
#   Origional code available at: https://github.com/ScottFerson/pba.r
######


#ii()    = range(0, stop = (ProbabilityBoundsAnalysis.steps-1)/ProbabilityBoundsAnalysis.steps, length = ProbabilityBoundsAnalysis.steps)

ii()    = [        0; collect((1:(parametersPBA.steps-1)) / parametersPBA.steps)];
iii()   = [  parametersPBA.bOt; collect((1:(parametersPBA.steps-1)) / parametersPBA.steps)];
jj()    = [           collect((1:(parametersPBA.steps-1)) / parametersPBA.steps); 1 ];
jjj()   = [           collect((1:(parametersPBA.steps-1)) / parametersPBA.steps); parametersPBA.tOp ];

##
# Should be able to env an array of p-boxes.
##
"""
    env(x :: pbox, y :: pbox)

Envelope. Returns the union of pboxes.

# Examples
```jldoctest
julia> a = U(0,1);

julia> b = U(1,2);

julia> c = env(a, b)
Pbox: 	  ~ uniform ( range=[0.0, 2.0], mean=[0.5, 1.5], var=0.083333)
```
See also: [`imp`](@ref), [`makepbox`](@ref)
"""
function env(x :: UncertainNumber, y :: UncertainNumber)
    x = makepbox(x);
    y = makepbox(y);

    na = x.name;
    sh = x.shape;
    bounded = x.bounded;

    u = min.(x.u, y.u);
    d = max.(x.d, y.d);
    ml = min(x.ml, y.ml);
    mh = max(x.mh, y.mh);
    vl = min(x.vl, y.vl);
    vh = max(x.vh, y.vh);
    bounded[1] = min(x.bounded[1], y.bounded[1]);
    bounded[2] = min(x.bounded[2], y.bounded[2]);

    if (y.name != na) na = "";end
    if (y.shape != sh) sh = "";end

    return pbox(u, d, ml=ml, mh=mh, vl=vl, vh=vh, name=na, shape=sh, bounded = bounded)

end
env(a...) = reduce(env, a)
env(a::Vector{UncertainNumber}) = reduce(env, a)

∪(x :: pbox, y :: UncertainNumber) = env(x, y)
∪(x :: UncertainNumber, y :: pbox) = env(x, y)

###
#   Computes the imprint (imposition)
###
"""
    imp(x :: pbox, y :: pbox)

Imprint. Returns the intersection of pboxes. Any number of pboxes may be input

# Examples
```jldoctest
julia> a = U(interval(0, 1), 2);

julia> b = U(1, 2);

julia> c = imp(a, b)
Pbox: 	  ~  ( range=[1.0, 2.0], mean=1.5, var=0.083333)
```
See also: [`env`](@ref), [`makepbox`](@ref)
"""
function imp(x...; naRm = false )
    elts = makepbox(x...);
    u = elts[1].u;
    d = elts[1].d;
    ml = elts[1].ml;
    mh = elts[1].mh;
    vl = elts[1].vl;
    vh = elts[1].vh;
    na = elts[1].name;
    sh = elts[1].shape;
    bounded =elts[1].bounded;

    for i =2:length(elts)
        u = max.(u,elts[i].u);
        d = min.(d,elts[i].d);
        ml = max(ml,elts[i].ml);
        mh = min(mh,elts[i].mh);
        vl = max(vl,elts[i].vl);
        vh = min(vh,elts[i].vh);
        bounded[1] = max(bounded[1], elts[i].bounded[1]);
        bounded[2] = max(bounded[2], elts[i].bounded[2]);
    end

    if (any(d[:] < u[:])) throw(ArgumentError("Imprint does not exist"));
    else
        return pbox(u, d, ml=ml, mh=mh, vl=vl, vh=vh, bounded = bounded)
    end

end
∩(x :: pbox, y :: UncertainNumber) = imp(x, y)
∩(x :: UncertainNumber, y :: pbox) = imp(x, y)


###
#   Constructs a pbox from a parametric distribution. Takes endpoints of given intervals.
#   Requires a dname function, which will give quantiles to the pbox.
#   envConstFunc is a version of this function which does not require a dname function
###
function envConstruct(dname, i , j, x...)

    a = env(
    map.(dname, left(i), left(j), x...),
    map.(dname, right(i), right(j), x...),
    map.(dname, left(i), right(j), x...),
    map.(dname, right(i), left(j), x...));

    return a;

end

###
#   Normal Distribution
###
function Snormal(mean = missing, std=missing; median=missing, mode=missing,
    cv=missing, iqr=missing, var=missing, name="", x...)

    if (ismissing(mean) && !ismissing(median))
        mean = median;
    end
    if (ismissing(mean) && !ismissing(mode))
        mean = mode;
    end
    if (ismissing(std) && !ismissing(cv) && !ismissing(mean))
        std = mean * cv;
    end
    if (!ismissing(mean) && !ismissing(std))
        return pbox(quantile.(Normal(mean,std),iii()), quantile.(Normal(mean,std),jjj()),
        shape="normal", name=name, ml=mean, mh=mean, vl=std^2, vh=std^2)
    else
        throw(ArgumentError("not enough information to specify the normal distribution"));
    end
end

"""
    normal(mean :: Interval, std :: Interval)

Normal shaped pbox. Parameters can be Real or Intervals.

# Constructors
* `normal`
* `N`
* `gaussian`

# Examples
```jldoctest
julia> a = normal(interval(0, 1), interval(1,2))
Pbox: 	  ~ normal ( range=[-6.1805, 7.1805], mean=[0.0, 1.0], var=[1.0, 4.0])
```
See also: [`uniform`](@ref), [`lognormal`](@ref), [`meanMinMax`](@ref), [`plot`](@ref)
"""
normal(mean = 0, std = 1, x...) = envConstruct(Snormal, mean, std, x...);
N = gaussian = normal

#Normal(mean :: Union{AbstractInterval,AbstractPbox}, std :: Union{AbstractInterval,AbstractPbox}, x...) = normal(mean, std, x...);

#Normal(mean, std) = normal(mean, std)
#Normal(mean, std, x...) = normal(mean, std, x...);

###
#   Uniform Distribtion
###


###
#   Env constructor specifically for the uniform distribtion.
###
function envUnif( i, j, x...)

    i, j, _, _ = checkMomentsAndRanges(i,j)
    a = env(
    map.(Suniform, left(i), left(j), 1, x...),
    map.(Suniform, right(i), right(j),x...),
    map.(Suniform, left(i), right(j), 1, x...),
    map.(Suniform, right(i), left(j),x...));
    return a
end

function Suniform(Min, Max, case = 1; name="")

    if (case==1) && (Max <= Min); Min = Max; end
    if Max <= Min; Max = Min; end

    m = (Min+Max)/2;
    v = (Min-Max)^2/12;
    if Max == Min; return pbox(Min,shape="uniform", name=name, ml=m, mh=m, vl=v, vh=v,bounded = [true, true]);end

    return (pbox(quantile.(Uniform(Min,Max),ii()), quantile.(Uniform(Min,Max),jj()),
    shape="uniform", name=name, ml=m, mh=m, vl=v, vh=v, bounded = [true, true]));

end

"""
    uniform(min :: Interval, max :: Interval)

Uniform shaped pbox. Parameters can be Real or Intervals.

# Constructors
* `uniform`
* `U`

# Examples
```jldoctest
julia> a = uniform(interval(0, 1), interval(1,2))
Pbox: 	  ~ uniform ( range=[0.0, 2.0], mean=[0.5, 1.5], var=[0.0, 0.33333])
```
See also: [`normal`](@ref), [`beta`](@ref), [`meanMinMax`](@ref), [`plot`](@ref)
"""
U = uniform(Min, Max, x...) = envUnif(Min,Max, x...)

function uniform(;mean, std)

    meanLo = left(mean); meanHi = right(mean)

    m = meanLo

    Min = m - sqrt(12)/2 * std
    Max = m + sqrt(12)/2 * std

    u1 = U(Min.hi, Max.lo)
    u2 = U(Min.lo, Max.hi)

    uu = env(u1, u2)

    m = meanHi

    Min = m - sqrt(12)/2 * std
    Max = m + sqrt(12)/2 * std

    u1 = U(Min.hi, Max.lo)
    u2 = U(Min.lo, Max.hi)

    uu = env(u1, u2, uu)

    return uu
end


###
#   Other Parametric distribution constructors
###

###
#   Constructs a pbox from a parametric distribution. Takes endpoints of given intervals.
#   Will work for all parametric functions up to 2 parameter. And if pbox is envelope of endpoints
#   of parameters
###
function envConstFunc(Dist, i, j, name, shape, Bounded)

    a = env(
    Sdist(Dist, left(i),  left(j),  name, shape, Bounded),
    Sdist(Dist, right(i), right(j), name, shape, Bounded),
    Sdist(Dist, left(i),  right(j), name, shape, Bounded),
    Sdist(Dist, right(i), left(j),  name, shape, Bounded));

    a.bounded = Bounded;
    return a;

end

function envConstFunc1(Dist, i, name, shape, Bounded)

    a = env(
    Sdist1(Dist, left(i),  name, shape, Bounded),
    Sdist1(Dist, right(i), name, shape, Bounded))

    return a;

end

###
#   The dname function for parametric distributions
###
function Sdist(DistFunc, i, j, name, shape, Bounded)

    is = iii(); js = jjj();
    if Bounded[1]; is = ii(); end
    if Bounded[2]; js = jj(); end

    Dist = DistFunc(i,j)
    m = mean(Dist); v = var(Dist);
    return pbox(quantile.(Dist,is), quantile.(Dist,js),shape=shape, name=name, ml=m, mh=m, vl=v, vh=v, bounded = Bounded)

end

function Sdist1(DistFunc, i, name, shape, Bounded)

    is = iii(); js = jjj();
    if Bounded[1]; is = ii(); end
    if Bounded[2]; js = jj(); end

    Dist = DistFunc(i)
    m = mean(Dist); v = var(Dist);
    return pbox(quantile.(Dist,is), quantile.(Dist,js),shape=shape, name=name, ml=m, mh=m, vl=v, vh=v, bounded = Bounded)

end


function Sdist3(DistFunc, i, j, k, name, shape, Bounded)

    is = iii(); js = jjj();
    if Bounded[1]; is = ii(); end
    if Bounded[2]; js = jj(); end

    Dist = DistFunc(i, j, k)
    m = mean(Dist); v = var(Dist);
    return pbox(quantile.(Dist,is), quantile.(Dist,js),shape=shape, name=name, ml=m, mh=m, vl=v, vh=v, bounded = Bounded)

end


###
#   For informaton about what parameters are, see the Distributions.jl package
###

"""
    beta(α :: Interval, β :: Interval)

Beta shaped pbox. Parameters can be Real or Intervals.

# Examples
```jldoctest
julia> a = beta(2,interval(3,4))
Pbox: 	  ~ beta ( range=[0.0, 1.0], mean=[0.33333, 0.4], var=[0.031746, 0.04])
```
See also: [`KN`](@ref), [`meanMinMax`](@ref), [`plot`](@ref)
"""
function beta( α=1,    β=1,    name = "")
    if α == 0 && β == 0; return pbox(0,1);end
    if α == 0; return pbox(0); end
    if β == 0; return pbox(1);end
    return envConstFunc( Beta, α, β, name, "beta", [true, true])
end

betaPrime(  α=1,    β=1,    name = "")      = envConstFunc(     BetaPrime,  α, β, name, "betaPrime",    [true,false])   # http://en.wikipedia.org/wiki/Beta_prime_distribution
biweight(   α,      β,      name = "")      = envConstFunc(     Biweight,   α, β, name, "biweight",     [true, true])
cauchy(     α=0,    β=1,    name = "")      = envConstFunc(     Cauchy,     α, β, name, "cauchy",       [false, false]) # http://en.wikipedia.org/wiki/Cauchy_distribution
chi(        k,              name = "")      = envConstFunc1(    Chi,        k,    name, "chi",          [true,false])   # http://en.wikipedia.org/wiki/Chi_distribution
chisq(      k,              name = "")      = envConstFunc1(    Chisq,      k,    name, "chusq",        [true,false])   # http://en.wikipedia.org/wiki/Chi-squared_distribution
cosine(     α,      β,      name = "")      = envConstFunc(     Cosine,     α, β, name, "cosine",       [true,true])    # http://en.wikipedia.org/wiki/Raised_cosine_distribution
epanechnikov(α,     β,      name = "")      = envConstFunc(     Epanechnikov, α, β, name, "epanechnikov",[true,true])
exponential(θ=1,            name = "")      = envConstFunc1(    Exponential, θ,   name, "exponential",  [true,false])   # http://en.wikipedia.org/wiki/Exponential_distribution
fDist(      ν1,     ν2,     name = "")      = envConstFunc(     FDist,     ν1,ν2, name, "fDist",        [true,false])   # http://en.wikipedia.org/wiki/F-distribution
gamma(      α=1,    θ=1,    name = "")      = envConstFunc(     Gamma,      α, θ, name, "gamma",        [true,false])   # http://en.wikipedia.org/wiki/Gamma_distribution
laplace(    μ=0,    θ=1,    name = "")      = envConstFunc(     Laplace,    μ, θ, name, "laplace",      [false,false])  # http://en.wikipedia.org/wiki/Laplace_distribution
levy(       μ=0,    θ=1,    name = "")      = envConstFunc(     Levy,       μ, θ, name, "levy",         [false,false])  # http://en.wikipedia.org/wiki/Laplace_distribution

function pbaLogNormal(m, std)
    γ = 1+std^2/m^2
    μ = log(m/sqrt(γ))
    σ = sqrt(log(γ))
    return LogNormal(μ,σ)
end

"""
    lognormal(μ :: Interval, std :: Interval)

Lognormal shaped pbox. Parameters can be Real or Intervals.

See also: [`KN`](@ref), [`meanMinMax`](@ref)
"""
lognormal(  μ=1,    θ=1,    name = "")      = envConstFunc(     pbaLogNormal,  μ, θ, name, "lognormal",    [true,false])   # http://en.wikipedia.org/wiki/Log-normal_distribution

###############################
#   Confidence boxes
###############################

"""
    KN(k :: Interval, n :: Interval)

k out of N confidence box (c-box), a pbox shaped confidence structure. Quantifies inferential uncertainty in binomial counts,
where k successes were observed out of n trails. One sided or two sided confidence intervals may be drawn

# Constructors
* `KN`
* `kn`

See also: [`meanMinMax`](@ref), [`plot`](@ref)
"""
function KN(k, n)
    if left(k) < 0
        throw(ArgumentError("k must be greater than 0, provided k = $k"))
    end
    if right(n) < right(k)
        throw(ArgumentError("k < n must be true, provided: k = $k | n = $n"))
    end
    return env( beta( left(k), right(n)-left(k)+1 ),  beta( right(k)+1, max(0, left(n)-right(k)) ) )
end
kn(x...) = KN(x...)

###############################
#   Distribution-free constructors for pboxes
###############################


###
#   Chebyshev inequality and mean - var pbox
###

function chebIneq(mean = 0, var = 1, name = "")
    p = iii();
    u = mean .- sqrt(var) .* sqrt.(1 ./p .-1)
    p = jjj();
    d = mean .+ sqrt(var) .* sqrt.(p ./(1 .-p))
    return pbox(u, d, shape = "chebyshev", name = name, ml = mean, mh = mean, vl = var, vh = var)
end

function meanVar(mean = 0, var = 1, name = "")

    if !isa(mean, Interval) && !isa(var, Interval);
        return chebIneq(mean, var,name)
    end
     Envelope = env(
        chebIneq(left(mean), left(var)),
        chebIneq(left(mean), right(var)),
        chebIneq(right(mean), left(var)),
        chebIneq(right(mean), right(var)),
    )
    Envelope.shape = "chebyshev"; Envelope.name = name;
    return Envelope
end

Chebyshev(x...) = meanVar(x...);    chebyshev(x...) = meanVar(x...);
cheb(x...)      = meanVar(x...);    MeanVar(x...)   = meanVar(x...);
meanvar(x...)    = meanVar(x...);

meanStd(mean,std, name = "") = meanVar(mean, std^2, name)

###
#   Markov inequality and mean - min pbox
###
function markovIneq(mean = 1, min = 0, name = "")

    if mean < min; min = mean; end
    p = jjj(); numS = length(p);
    d = ( (mean - min) ./(1 .- p) ) .+ min;
    return pbox(ones(numS)*min, d, shape = "markov", name = name, ml = mean, mh = mean, vl =0 , vh = Inf, bounded = [true, false])
end

function meanMin(mean = 1, min = 0, name = "")

    if right(mean) < left(min); throw(ArgumentError("Inconsistent arguements: \n            $mean < $min"));end

    if !isa(mean, Interval) return markovIneq(mean, left(min), name); end

    Envelope = env(
        markovIneq(left(mean), left(min)),
        markovIneq(right(mean), left(min)),
    )
    Envelope.shape = "markov"; Envelope.name = name;
    return Envelope
end

Markov(x...)    = meanMin(x...);    markov(x...)    = meanMin(x...);
MeanMin(x...)   = meanMin(x...);    meanmin(x...)   = meanMin(x...);

meanMax(mean = 0, max = 1, name = "") = pbox(negate(meanMin(- (mean), - (max))) , name = name)
MeanMax(x...) = meanMax(x...);      meanmax(x...)= meanMax(x...);

###
#   Cantelli inequality and mean - min - max pbox
###
function cantelliIneq(mean = 0.5, min = 0, max = 1, name = "")

    if mean < min; min = mean; end
    if mean > max; max = mean; end

    mid = (max - mean) / (max - min);
    p = ii();
    u = ifelse.(p .<= mid, min, (mean - max) ./p .+ max )
    p = jj();
    d = ifelse.(mid .<= p, max, (mean .- min .*p ) ./(1 .- p))

    return pbox(u, d, shape = "cantelli", name = name, ml = mean, mh = mean, vl = 0,
    vh = (max-min)*(max-mean)-(max-mean)*(max-mean), bounded = [true, true])

end

function meanMinMax(mean = 0.5, min = 0, max = 1, name = "")

    min, max, mean1, var = checkMomentsAndRanges(min,max,mean);

    if !isa(mean, Interval) return cantelliIneq(mean, left(min), right(max), name); end

    Envelope = env(
        cantelliIneq(left(mean1), left(min), right(max)),
        cantelliIneq(right(mean1), left(min), right(max))
    )
    Envelope.shape = "cantelli"; Envelope.name = name;
    Envelope.vl = var.lo; Envelope.vh = var.hi;
    return Envelope
end

Cantelli(x...)      = meanMinMax(x...);     cantelli(x...)      = meanMinMax(x...);
MeanMinMax(x...)    = meanMinMax(x...);     meanminmax(x...)    = meanMinMax(x...);


function cantelliIneq2(minX = 0, meanX = 0.5, varX = 1, name = "")

    if meanX < minX; minX = meanX; end

    meanY = meanX - minX; # Transform to 0

    p = ii();
    u = max.(0, meanY .-  sqrt.(varX .*(1 ./p .- 1)))
    p = jjj()
    d = Base.min.(meanY ./ (1 .- p), meanY .+ sqrt.((p .* varX) ./ (1 .- p)))

    s = pbox(u, d, shape = "cantelli2", name = name, ml = meanY, mh = meanY, vl = varX,
    vh = varX, bounded = [true, false])

    return s + minX
end

function minMeanVar(min = 0, mean = 0.5, var = 1, name = "")

    mean = interval(mean); var = interval(var);

    Envelope = env(
        cantelliIneq2(left(min), left(mean), right(var)),
        cantelliIneq2(right(min), right(mean), right(var))
    )
    Envelope.shape = "cantelli"; Envelope.name = name;
    Envelope.vl = var.lo; Envelope.vh = var.hi;
    return Envelope
end

function maxMeanVar(max = 1, mean = 0.5, var=1, name = "")

    mean = interval(mean); var = interval(var);

    minZ = -max
    meanZ = - mean

    Envelope = env(
        cantelliIneq2(left(minZ), left(meanZ), right(var)),
        cantelliIneq2(right(minZ), right(meanZ), right(var))
    )
    Envelope.shape = "cantelli"; Envelope.name = name;
    Envelope.vl = var.lo; Envelope.vh = var.hi;
    return -Envelope

end

###
#   Ferson Pbox (min - max - mean - var)
###
function mmms(minimum :: Real, maximum :: Real, mean :: Real, stddev :: Real; steps = 200)

    if iszero(maximum - minimum); return pbox((maximum + minimum)/2); end
    range = maximum - minimum;
    min1, max1, mean1,var1 = checkMomentsAndRanges(minimum,maximum, mean, stddev^2);

    s = sqrt(var1)

    ml = (left(mean1)-minimum)/range;       sl = left(s)/range;
    mr = (right(mean1)-minimum)/range;      sr = right(s)/range;

    n  = steps;
    us = zeros(n);
    ds = zeros(n);
    for i = 1:n

        p = (i-1)/n
        if p <= 0.0; x2 = 0.0;
        else; x2 = ml - sr * sqrt(1 / p - 1); end
        if ml + p <= 1.0;
            x3 = 0.0;
        else
            x5 = p*p + sl*sl - p;
            if x5 >= 0.0
                x4 = 1.0 - p + sqrt(x5);
                if x4 < ml
                    x4 = ml;
                end
            else
                x4 = ml;
            end
            x3 = (p + sl*sl + x4*x4 - 1.0) / (x4 + p - 1.0);
        end

        if ((p <= 0.0) || (p <= (1.0 - ml)))
            x6 = 0.0;
        else
            x6 = (ml - 1.0) / p + 1.0;
        end

        us[i] = max(max(max(x2, x3), x6),0.0) * range + minimum;

        p = i/n;

        if (p >= 1.0)
            x2 = 1.0;
        else
            x2 = mr + sr * sqrt(1.0/(1.0/p - 1.0));
        end
        if (mr + p >= 1.0)
            x3 = 1.0;
        else
            x5 = p*p + sl*sl - p;
            if x5 >= 0.0
                x4 = 1.0 - p - sqrt(x5);
                if (x4 > mr) x4 = mr; end

            else
                x4 = mr;
            end
            x3 = (p + sl*sl + x4*x4 - 1.0) / (x4 + p - 1.0) - 1.0;
        end
        if (((1.0 - mr) <= p) || (1.0 <= p))
            x6 = 1.0;
        else
            x6 = mr / (1.0 - p);
        end

        ds[i] = min(min(min(x2,x3),x6),1.0) * range + minimum;

    end

    return pbox(us, ds, ml = left(mean1), mh = right(mean1), vl = left(var1), vh = right(var1), shape = "ferson")
end

#mmms(minimum :: Real, maximum :: Real, mean :: Real, std :: Real; steps = 200) = mmms(minimum, maximum, interval(mean), interval(std), steps = steps)

mmmv(minimum :: Real, maximum :: Real, mean :: Real, var :: Real; steps = 200) = mmms(minimum, maximum, mean, sqrt(var), steps = steps)
minMaxMeanVar(minimum :: Real, maximum :: Real, mean :: Real, var :: Real; steps = 200) = mmmv(minimum , maximum , mean, var; steps = steps)
minMaxMeanStd(minimum :: Real, maximum :: Real, mean :: Real, std :: Real; steps = 200) = mmmv(minimum , maximum , mean, std^2; steps = steps)


Ferson(x...)        = minMaxMeanVar(x...);         ferson(x...) = minMaxMeanVar(x...);
MinMaxMeanVar(x...) = minMaxMeanVar(x...);  minmaxmeanvar(x...) = minMaxMeanVar(x...);


# This cut could be quicker, and should also allow interval arguments
"""
    cut(x :: pbox, p :: Real)

returns a vertical cut of a pbox at cdf value p, for p ∈ [0, 1]

# Constructors
* `cut(x :: pbox, p :: Real)`
* `cut(x :: pbox, p :: Interval)`

See also: [`rand`](@ref), [`cdf`](@ref), [`mass`](@ref)
"""
function cut(x, p :: Real; tight :: Bool = true)

    x = makepbox(x);

    if (p<0 || p>1) throw(ArgumentError("Second argument must be a probability between zero and one")); end

    if (p < parametersPBA.bOt && !x.bounded[1]); return interval(-∞, x.d[1]);end
    if (p > parametersPBA.tOp && !x.bounded[2]); return interval(x.u[end], ∞);end

    long = x.n;
    if (tight) return (interval(x.u[Int(min.(long,(mod(p*long, 1)==0)+ceil(p*x.n)))], x.d[ Int(max.(1,ceil(p*x.n)))])) end
    if (p == 1) lower = long; elseif (mod(p,(1/long)) == 0) lower = round(p*long); else lower = ceil(p*long); end
    if (p == 0) upper = 1; elseif (mod(p,(1/long)) == 0) upper = round(p*long)+1; else upper = floor(p*long)+1; end
    return interval(x.u[Int(max(lower,1))], x.d[Int(min(upper,long))]);

end

cut(x, p :: Interval; tight :: Bool = true) = interval(cut(x,left(p)), cut(x,right(p)))




"""
    rand(x :: pbox, n :: Int64)

returns n number of random intervals from pbox x

# Examples
```jldoctest
julia> a = N(0,1);

julia> rand(a,5)
5-element Vector{Interval{Float64}}:
  [0.538836, 0.553385]
 [-0.331854, -0.318639]
  [0.125661, 0.138305]
 [-1.40508, -1.3722]
 [-0.31864, -0.30548]
```
See also: [`cut`](@ref), [`cdf`](@ref), [`mass`](@ref)
"""
rand(a :: pbox, n :: Int64 = 1) = cut.(a,rand(n));


###
#   Checks consitency of provided moments and ranges, which may all be intervals.
#   Returns consitent intervals. If var not provided, returns best var. Could maybe make this a macro?
###
function checkMomentsAndRanges(Min, Max, Mean=interval(left(Min),right(Max)), Var = interval(0,Inf))

    if right(Max) < left(Min); throw(ArgumentError("Inconsistent bounds max < min: $Max < $Min")); end

    range = interval(left(Min), right(Max));
    if Mean ∩ range == ∅;
        throw(ArgumentError("Inconsistent arguements: mean ∩ [min, max] = ∅ \n            $Mean ∩ $range = ∅"))
    end

    # Fix min and max
    MAXl = max(left(Max), left(Min), left(Mean))
    MINh = min(right(Min), right(Max), right(Mean))

    Max = interval(MAXl,  right(Max));
    Min = interval(left(Min), MINh)

    # Fix mean
    mh = min(right(Mean), right(Max));
    ml = max(left(Mean), left(Min));

    Mean = interval(ml,mh);

    # Find variance range
    MIN = left(Min); MAX = right(Max);
    v1 = (MAX-MIN)*(MAX-ml)-(MAX-ml)*(MAX-ml);
    v2 = (MAX-MIN)*(MAX-mh)-(MAX-mh)*(MAX-mh);
    v3 = 0;

    # Use mid point if it is in mean
    mid = (MAX-MIN)/2;
    if (mid ∈ Mean); v3 = (MAX-MIN)*(MAX-mid)-(MAX-mid)*(MAX-mid); end

    vh = max(v1, v2, v3); vl = 0;
    maxVar = interval(vl,vh);

    if (Var ∩ maxVar == ∅); throw(ArgumentError("Provided information not valid. Variance ∩ VarBounds = ∅.\n       $Var ∩ $maxVar = ∅")); end

    if !(Var ⊆ maxVar)  # If it's already a subset no need to bother
        vl = max(left(Var), left(maxVar));
        vh = min(right(Var), right(maxVar));
        Var = interval(vl, vh)
    end

    return Min, Max, Mean, Var

end
