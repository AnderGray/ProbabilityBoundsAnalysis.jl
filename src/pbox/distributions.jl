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

function envConst(dname, i , j,x...)

    a = env(
    map.(dname, left(i), left(j),x...),
    map.(dname, right(i), right(j),x...),
    map.(dname, left(i), right(j),x...),
    map.(dname, right(i), left(j),x...));
    
    a.dids = "PB $(uniquePbox())";
    return a;

end

###
#   Normal Distribution
###

function Snormal0(normmean, normstd, name="")
    return pbox(quantile.(Normal(normmean,normstd),iii()), quantile.(Normal(normmean,normstd),jjj()),
    shape="normal", name=name, ml=normmean, mh=normmean, vl=normstd^2, vh=normstd^2)
end

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
        return Snormal0(mean,std,name)
    else
        throw(ArgumentError("not enough information to specify the normal distribution"));
    end
end

N = normal = gaussian(mean, std, x...) = envConst(Snormal, mean, std, x...);
Normal(mean :: Union{AbstractInterval,AbstractPbox}, std :: Union{AbstractInterval,AbstractPbox}, x...) = normal(mean, std, x...);

#Normal(mean, std) = normal(mean, std)
#Normal(mean, std, x...) = normal(mean, std, x...);

###
#   Uniform Distribtion
###
function Suniform(min, max, case = 1, name="")

    if (case==1) && (max <= min); min = max; end
    if max <= min; max = min; end

    m = (min+max)/2;
    v = (min-max)^2/12;
    if max == min; return pbox(min,shape="uniform", name=name, ml=m, mh=m, vl=v, vh=v);end

    return (pbox(quantile.(Uniform(min,max),ii()), quantile.(Uniform(min,max),jj()),
    shape="uniform", name=name, ml=m, mh=m, vl=v, vh=v));

end


function envUnif(i , j,x...)

    a = env(
    map.(Suniform, left(i), left(j), 1, x...),
    map.(Suniform, right(i), right(j),x...),
    map.(Suniform, left(i), right(j), 1, x...),
    map.(Suniform, right(i), left(j),x...));

end

U = uniform(min, max, x...) = envUnif(min,max, x...)



###
#   Uniform Distribtion
###

function Sbeta(α, β,name="")

    a = Beta(α, β);
    m = mean(a); v = var(a);

    return pbox(quantile.(a,ii()), quantile.(a,jj()),shape="beta", name=name, ml=m, mh=m, vl=v, vh=v)

end

beta(α, β, x...) = envConst(Sbeta, α, β)

###
#   Distribution-free constructors for pboxes
###



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
