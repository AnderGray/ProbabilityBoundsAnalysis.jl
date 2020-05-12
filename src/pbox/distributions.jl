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
#   Port of R code pba.r by Scott Ferson and Jason O'Rawe, Applied Biomathematics
#   Origional code available at: https://github.com/ScottFerson/pba.r
######


#ii()    = range(0, stop = (ProbabilityBoundsAnalysis.steps-1)/ProbabilityBoundsAnalysis.steps, length = ProbabilityBoundsAnalysis.steps)

ii()    = [        0; collect((1:(ProbabilityBoundsAnalysis.steps-1)) / ProbabilityBoundsAnalysis.steps)];
iii()   = [  ProbabilityBoundsAnalysis.bOt; collect((1:(ProbabilityBoundsAnalysis.steps-1)) / ProbabilityBoundsAnalysis.steps)];
jj()    = [           collect((1:(ProbabilityBoundsAnalysis.steps-1)) / ProbabilityBoundsAnalysis.steps); 1 ];
jjj()   = [           collect((1:(ProbabilityBoundsAnalysis.steps-1)) / ProbabilityBoundsAnalysis.steps); ProbabilityBoundsAnalysis.tOp ];

##
# Should be able to env an array of p-boxes. Feature to be added 
##
function env(x...; naRm = false )
    elts = makepbox(x...);
    u = elts[1].u;
    d = elts[1].d;
    ml = elts[1].ml;
    mh = elts[1].mh;
    vl = elts[1].vl;
    vh = elts[1].vh;
    dids = elts[1].dids;
    na = elts[1].name;
    sh = elts[1].shape;

    for i =2:length(elts)
        u = min.(u,elts[i].u);
        d = max.(d,elts[i].d);
        ml = min(ml,elts[i].ml);
        mh = max(mh,elts[i].ml);
        vl = min(vl,elts[i].vl);
        vh = max(vh,elts[i].vh);
        dids = "$dids $(elts[i].dids)"
        if (elts[i].name != na) na = "";end
        if (elts[i].shape != sh) sh = "";end
    end
    return pbox(u, d, ml=ml, mh=mh, vl=vl, vh=vh, dids=dids, name=na, shape=sh)

end


###
#   Computes the imprint (imposition)
###
function imp(x...; naRm = false )
    elts = makepbox(x...);
    u = elts[1].u;
    d = elts[1].d;
    ml = elts[1].ml;
    mh = elts[1].mh;
    vl = elts[1].vl;
    vh = elts[1].vh;
    dids = elts[1].dids;
    na = elts[1].name;
    sh = elts[1].shape;

    for i =2:length(elts)
        u = max.(u,elts[i].u);
        d = min.(d,elts[i].d);
        ml = max(ml,elts[i].ml);
        mh = min(mh,elts[i].mh);
        vl = max(vl,elts[i].vl);
        vh = min(vh,elts[i].vh);
        dids = "$dids $(elts[i].dids)"
    end

    if (any(d[:] < u[:])) throw(ArgumentError("Imprint does not exist"));
    else
        return pbox(u, d, ml=ml, mh=mh, vl=vl, vh=vh, dids=dids)
    end

end


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
    
    a.dids = "PB $(uniquePbox())";
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

    a = env(
    map.(Suniform, left(i), left(j), 1, x...),
    map.(Suniform, right(i), right(j),x...),
    map.(Suniform, left(i), right(j), 1, x...),
    map.(Suniform, right(i), left(j),x...));

end

function Suniform(min, max, case = 1; name="")

    if (case==1) && (max <= min); min = max; end
    if max <= min; max = min; end

    m = (min+max)/2;
    v = (min-max)^2/12;
    if max == min; return pbox(min,shape="uniform", name=name, ml=m, mh=m, vl=v, vh=v);end

    return (pbox(quantile.(Uniform(min,max),ii()), quantile.(Uniform(min,max),jj()),
    shape="uniform", name=name, ml=m, mh=m, vl=v, vh=v));

end

U = uniform(min, max, x...) = envUnif(min,max, x...)


###
#   Other Parametric distribution constructors
###

#=
function Sbeta(α, β, name="")

    a = Beta(α, β);
    m = mean(a); v = var(a);

    return pbox(quantile.(a,ii()), quantile.(a,jj()), shape="beta", name=name, ml=m, mh=m, vl=v, vh=v)

end

beta(α, β, x...) = envConst(Sbeta, α, β, x...)
=#


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
    
    a.dids = "PB $(uniquePbox())";
    return a;

end

function envConstFunc1(Dist, i, name, shape, Bounded)

    a = env(
    Sdist1(Dist, left(i),  name, shape, Bounded),
    Sdist1(Dist, right(i), name, shape, Bounded))
    
    a.dids = "PB $(uniquePbox())";
    return a;

end

function envConstFunc3(Dist, i, j, k,  name, shape, Bounded)

    a = env(
    Sdist(Dist, left(i),  left(j),  name, shape, Bounded),
    Sdist(Dist, right(i), right(j), name, shape, Bounded),
    Sdist(Dist, left(i),  right(j), name, shape, Bounded),
    Sdist(Dist, right(i), left(j),  name, shape, Bounded));
    
    a.dids = "PB $(uniquePbox())";
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
    return pbox(quantile.(Dist,is), quantile.(Dist,js),shape=shape, name=name, ml=m, mh=m, vl=v, vh=v)

end

function Sdist1(DistFunc, i, name, shape, Bounded)

    is = iii(); js = jjj(); 
    if Bounded[1]; is = ii(); end
    if Bounded[2]; js = jj(); end

    Dist = DistFunc(i)
    m = mean(Dist); v = var(Dist);
    return pbox(quantile.(Dist,is), quantile.(Dist,js),shape=shape, name=name, ml=m, mh=m, vl=v, vh=v)

end


function Sdist3(DistFunc, i, j, k, name, shape, Bounded)

    is = iii(); js = jjj(); 
    if Bounded[1]; is = ii(); end
    if Bounded[2]; js = jj(); end

    Dist = DistFunc(i, j, k)
    m = mean(Dist); v = var(Dist);
    return pbox(quantile.(Dist,is), quantile.(Dist,js),shape=shape, name=name, ml=m, mh=m, vl=v, vh=v)

end


###
#   For informaton about what parameters are, see the Distributions.jl package
###

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
erlang(     α=1,    β=1,    name = "")      = envConstFunc(     Erlang,     α, β, name, "erlang",       [true,false])   # http://en.wikipedia.org/wiki/Erlang_distribution
exponential(θ=1,            name = "")      = envConstFunc1(    Exponential, θ,   name, "exponential",  [true,false])   # http://en.wikipedia.org/wiki/Exponential_distribution         
fDist(      ν1,     ν2,     name = "")      = envConstFunc(     FDist,     ν1,ν2, name, "fDist",        [true,false])   # http://en.wikipedia.org/wiki/F-distribution
frechet(    α=1,    θ=1,    name = "")      = envConstFunc(     Frechet,    α, θ, name, "frechet",      [true,false])   # http://en.wikipedia.org/wiki/Fréchet_distribution
gamma(      α=1,    θ=1,    name = "")      = envConstFunc(     Gamma,      α, θ, name, "gamma",        [true,false])   # http://en.wikipedia.org/wiki/Gamma_distribution
ksdist(     n,              name = "")      = envConstFunc1(    KSDist,     n,    name, "ksdist",       [true,true])
laplace(    μ=0,    θ=1,    name = "")      = envConstFunc(     Laplace,    μ, θ, name, "laplace",      [false,false])  # http://en.wikipedia.org/wiki/Laplace_distribution
levy(       μ=0,    θ=1,    name = "")      = envConstFunc(     Levy,       μ, θ, name, "levy",         [false,false])  # http://en.wikipedia.org/wiki/Laplace_distribution
lognormal(  μ=0,    θ=1,    name = "")      = envConstFunc(     lognormal,  μ, θ, name, "lognormal",    [true,false])   # http://en.wikipedia.org/wiki/Log-normal_distribution

###############################
#   Confidence boxes
###############################

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
minvar(x...)    = meanVar(x...);

meanStd(mean,std, name = "") = meanVar(mean, std^2, name)

###
#   Markov inequality and mean - min pbox
###
function markovIneq(mean = 1, min = 0, name = "")

    if mean < min; min = mean; end
    p = jjj(); numS = length(p);
    d = ( (mean - min) ./(1 .- p) ) .+ min;
    return pbox(ones(numS)*min, d, shape = "markov", name = name, ml = mean, mh = mean, vl =0 , vh = Inf)
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
    
    return pbox(u, d, shape = "cantelli", name = name, ml = mean, mh = mean, vl = 0, vh = (max-min)*(max-mean)-(max-mean)*(max-mean))

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


###
#   Ferson Pbox (min - max - mean - var)
###
function FersonEvalEasy(min = 0, max = 1,  mean = 0.5, var = 0.1, name = "")
    
    fer = pbox(imp(meanVar(mean,var), imp(meanmin(mean,min),meanmax(mean,max))),
    ml = left(mean), mh = right(mean), vl = left(var), vh =right(var))
    fer.shape = "ferson"; fer.name = name;

    return fer
end

function minMaxMeanVar(min = 0, max = 1, mean = 0.5, var = 0.1, name = "")

    min1,max1,mean1,var1 = checkMomentsAndRanges(min,max,mean,var);

    if all(!isa([mean, var], Interval));
        return FersonEvalEasy(left(min),right(max),mean,var,name);
    end

    Envelope = env(
        FersonEvalEasy(left(min1), right(max1), left(mean1), right(var1)),            # pba.r has it like this (left(var) not used?)
        FersonEvalEasy(left(min1), right(max1), right(mean1), right(var1)),
    )
    Envelope.shape = "ferson"; Envelope.name = name;
    return Envelope

end

Ferson(x...)        = minMaxMeanVar(x...);         ferson(x...) = minMaxMeanVar(x...);
MinMaxMeanVar(x...) = minMaxMeanVar(x...);  minmaxmeanvar(x...) = minMaxMeanVar(x...);


# This cut could be quicker, and should also allow interval arguments
function cut(x, p :: Real; tight :: Bool = true)

    x = makepbox(x);
    if (p<0 || p>1) throw(ArgumentError("Second argument must be a probability between zero and one")); end
    long = x.n;
    if (tight) return (interval(x.u[Int(min.(long,(mod(p*long, 1)==0)+ceil(p*x.n)))], x.d[ Int(max.(1,ceil(p*x.n)))])) end
    if (p == 1) lower = long; elseif (mod(p,(1/long)) == 0) lower = round(p*long); else lower = ceil(p*long); end
    if (p == 0) upper = 1; elseif (mod(p,(1/long)) == 0) upper = round(p*long)+1; else upper = floor(p*long)+1; end
    return interval(x.u[Int(max(lower,1))], x.d[Int(min(upper,long))]);

end

cut(x, p :: Interval; tight :: Bool = true) = interval(cut(x,left(p),tight=tight), cut(x,right(p),tight=tight))

rand(a :: pbox, n :: Int64 = 1) = cut.(a,rand(n));


###
#   Checks consitency of provided moments and ranges, which may all be intervals. 
#   Returns consitent intervals. If var not provided, returns best var. Could maybe make this a macro?
###
function checkMomentsAndRanges(Min, Max, Mean=interval(left(min),right(max)), Var = interval(0,Inf))

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

#=
myvectormax <- function(..., na.rm = FALSE) {
  a <- c(...)
  if (((isTRUE(all.equal(0,length(which.max(a)))))) || (isTRUE(all.equal(0,length(a))))) -Inf else if (!na.rm && any(is.na(a))) NA else a[[which.max(a)]]
  }

myvectormin <- function(..., na.rm = FALSE) {
  a <- c(...)
  if (((isTRUE(all.equal(0,length(which.min(a)))))) || (isTRUE(all.equal(0,length(a))))) +Inf else if (!na.rm && any(is.na(a))) NA else a[[which.min(a)]]
  }
  =#
