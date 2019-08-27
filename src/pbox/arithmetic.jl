######
# This file is part of the pba.jl package.
#
# Definition of arithmetic between pboxes
#
#           University of Liverpool
######



###########################
# Convolutions operations #
###########################

function conv(x::Real, y::Real, op = +)

    x = makepbox(x);
    y = makepbox(y);

    if (op == -) return conv(x, negate(y), +);end
    if (op == /) return conv(x, reciprocate(y), *);end

    m = x.n;
    p = y.n;
    n = min(pba.steps, m*p);
    L = m * p / n;
    c = zeros(m*p);
    Zu = ones(n);
    Zd = ones(n);

    k = 1:n;

    Y = repeat(y.d[:], inner = m);
    for i=1:p
        c[m*(i-1)+1:m*i] = map(op, x.d[:], Y[m*(i-1)+1:m*i]);
    end

    c = sort(c);
    Zd = c[Int.(((k .- 1) .* L) .+ L)];

    ll = c;

    Y = repeat(y.u[:], inner = m);
    c = zeros(length(Y));
    for i=1:p
        c[m*(i-1)+1:m*i] = map(op, x.u[:], Y[m*(i-1)+1:m*i]);
    end

    c = sort(c);
    Zu = c[Int.(((k .- 1) .* L) .+ 1)];

    #println(Zu)
    #mean
    ml = -Inf;
    mh = Inf;

    if (op)∈([+,-,*])
        ml = map(op,x.ml,y.ml);
        mh = map(op,x.mh,y.mh);
    end

    # Variance
    vl = 0;
    vh = Inf

    if (op)∈([+,-])     # Where is *?
        vl = map(op,x.vl,y.vl);
        vh = map(op,x.vh,y.vh);
    end

    return pbox(Zu, Zd, ml = ml, mh = mh, vl=vl, vh=vh, dids="$(x.dids) $(y.dids)");
    #return ([Zu, Zd, ml, mh, vl, vh, "$(x.dids) $(y.dids)"]);

end

function convPerfect(x::Real, y::Real, op = +)
    if (op)∈([-,/])
        cu = map(op, x.u[:],y.d[:]);
        cd = map(op, x.d[:],y.u[:]);
    else
        cu = map(op, x.u[:], y.u[:]);
        cd = map(op, x.d[:], y.d[:]);
    end
    scu = sort(cu);
    scd = sort(cd);

    if (all(cu == scu) && all(cd == scd))
        return pbox(scu, scd,  dids="$(x.dids) $(y.dids) ", bob=x.bob)
    else return pbox(scu, scd,  dids="$(x.dids) $(y.dids) ")
      end
end

function convOpposite(x::Real, y::Real, op = +)

    if (op)∈([-,/])
        cu = map(op, x.u[:],y.d[end:-1:1]);
        cd = map(op, x.d[:],y.u[end:-1:1]);
    else
        cu = map(op, x.u[:], y.u[end:-1:1]);
        cd = map(op, x.d[:], y.d[end:-1:1]);
    end
    cu = sort(cu);
    cd = sort(cd);

    return pbox(cu, cd, dids = "$(x.dids) $(y.dids)")
end

function convFrechet(x::Real, y::Real, op = +)

    if (op == -) return (convFrechet(x,negate(y),+));end
    if (op == /) return (convFrechet(x,reciprocate(y),*));end
    #if (op==*) if (straddlingzero(x) || straddlingzero(y)) return (imp(balchprod(x,y),convFrechetNaive(x,y,*))); end
    ## Unsure about the above line. It looks like if it straddles 0, we need to do the naive frechet and the balch prod (?) and impose one on the other

    x = makepbox(x);
    y = makepbox(y);

    zu = zeros(pba.steps);
    zd = zeros(pba.steps);

    for i = 1:pba.steps

        j = i:pba.steps;
        k = pba.steps:-1:i;
        zd[i] = minimum(map(op, x.d[j],y.d[k]));

        j = 1:i;
        k = i:-1:1;
        zu[i] = maximum(map(op, x.u[j], y.u[k]));
    end
    #mean

    ml = -Inf;
    mh = Inf;
    if (op)∈([+,-])                 # We should be able to include * /  once we have momemnt prop
        ml = map(op,x.ml,y.ml)
        mh = map(op,x.mh,y.mh)
    end

    vl = 0;
    vh = Inf;

    if (op)∈([+,-])                 # Was commented below, can include onve mom prop is finished
        #zv = env(xv+yv-2* sqrt(xv*yv), xv+yv+2* sqrt(xv*yv))
        # vh <- x@v+y@v - 2*sqrt(x@v*y@v)
        # vl <- x@v+y@v + 2*sqrt(x@v*y@v)
    end

    return pbox(zu, zd, ml = ml, mh = mh, vl = vl, vh = vh, dids = "$(x.dids) $(y.dids) ");
end


function negate(x::pbox)
    if ((x.shape)∈(["uniform", "normal", "cauchy", "triangular"])) s = x.shape; else s = ""; end
    return pbox(-x.d[end:-1:1],-x.u[end:-1:1],shape=s,name = "", ml=-x.mh, mh=-x.ml, vl=x.vl, vh=x.vh, dids=x.dids, bob=oppositedep(x));
end

function complement(x::pbox)
    if ((x.shape)∈(["uniform", "normal", "cauchy", "triangular", "skew-normal"])) s = x.shape; else s = ""; end
    return pbox(1 .-x.d[end:-1:1],1 .-x.u[end:-1:1],shape=s,name = "", ml=1-x.mh, mh=1-x.ml, vl=x.vl, vh=x.vh, dids=x.dids, bob=oppositedep(x));
end


-(x::pbox) = negate(x);

+(x::AbstractPbox, y::AbstractPbox) = conv(x,y,+); # if(x==y) return 2*x; ????
-(x::AbstractPbox, y::AbstractPbox) = conv(x,y,-); # if(x==y) return 0;   ????
*(x::AbstractPbox, y::AbstractPbox) = conv(x,y,*); # if(x==y) return x^2; ????
/(x::AbstractPbox, y::AbstractPbox) = conv(x,y,/); # if(x==y) return 0;   ????

###
#   Conv of pboxes and intervals
###

# Probably will only need shift for + and - with reals
+(x :: AbstractPbox, y :: AbstractInterval) = conv(x,y,+);
+(x :: AbstractInterval, y :: AbstractPbox) = y + x;

-(x :: AbstractPbox, y :: AbstractInterval) = conv(x,y,-);
-(x :: AbstractInterval, y :: AbstractPbox) = -y + x;

*(x :: AbstractPbox, y :: AbstractInterval) = conv(x,y,*);
*(x :: AbstractInterval, y :: AbstractPbox) = y*x;

/(x :: AbstractPbox, y :: AbstractInterval) = conv(x,y,/);
/(x :: AbstractInterval, y :: AbstractPbox) = reciprocate(y) * x;

###
#   Conv of pboxes and reals
###

# Probably will only need shift for + and - with reals
+(x :: AbstractPbox, y :: Real) = conv(x,y,+);
+(x :: Real, y :: AbstractPbox) = y + x;

-(x :: AbstractPbox, y :: Real) = conv(x,y,-);
-(x :: Real, y :: AbstractPbox) = -y + x;

*(x :: AbstractPbox, y :: Real) = conv(x,y,*);
*(x :: Real, y :: AbstractPbox) = y*x;

/(x :: AbstractPbox, y :: Real) = conv(x,y,/);
/(x :: Real, y :: AbstractPbox) = reciprocate(y) * x;


oppositedep(x::pbox) = -x.bob;
perfectdep(x::pbox) = x.bob;
perfectopposite(m, x::pbox) = if (m<0) return oppositedep(x); else return perfectdep(x);end


###
#   Still needed
###

function reciprocate(x::pbox)

    if ((x.shape)∈(["Cauchy","{min, max, median}","{min, max, percentile}","{min, max}"]))  sh = x.shape;
    elseif (x.shape == "pareto")    sh = "power function";
    elseif (x.shape == "power function")    sh = "pareto";
    else sh = "";
    end

    return x;
end
function convFrechetNaive(x::pbox, y::pbox) return x;end
function imp(x::pbox, y:: pbox) return x; end
function balchProd(x::pbox, y::pbox) return x; end
