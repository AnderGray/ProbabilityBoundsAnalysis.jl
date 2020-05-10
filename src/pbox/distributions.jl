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

    if (any(d[:] < u[:])) throw(ArgumentError("Imposition does not exist"));
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

N = normal = gaussian(mean, std, x...) = envConstruct(Snormal, mean, std, x...);
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

###
#   Confidence boxes
###

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

###
#   Distribution-free constructors for pboxes
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


function meanMin(mean, min, name = "")
end

Markov(x...) = markov(x...) = MeanVar(x...) = meanVar(x...)

function meanMinMax(mean, min, max, name = "")
end

Cantelli(x...) = cantelli(x...) = MeanMinMax(x...) = meanminmax(x...) = meanMinMax(x...)

function meanMinMaxVar(mean, min, max, var, name = "")
end

Ferson(x...) = ferson(x...) = MeanMinMaxVar(x...) = meanminmaxvar(x...) = meanMinMaxVar(x...)


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
