######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Definition of arithmetic between pboxes
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
#
#   Port of R code pba.r by Scott Ferson and Jason O'Rawe, Applied Biomathematics
#   Origional code available at: https://github.com/ScottFerson/pba.r
#
######

#####
#   Known Bugs:
#
#   ->  Convolution with a copula defined (sigma convolution) not quite working as exected for negative correlations.
#		It could be the case that he imprints it with meanVarRange
#
#
####x




###########################
# Convolutions operations #
###########################

function conv(x::pbox, y::pbox; op = +, corr =0)
    if corr == 0; return convIndep(x,y,op = op); end
    if corr == 1; return convPerfect(x,y, op = op);end
    if corr ==-1; return convOpposite(x,y,op = op);end
    if corr ==interval(-1,1); return convFrechet(x,y, op = op); end

    if isinterval(corr)
        Lower = sigma(x, y, C = GauCopula(left(corr)), op=op);
        Upper = sigma(x, y, C = GauCopula(right(corr)), op=op);
        return env(Lower, Upper)
    end

    return sigma(x, y, C = GauCopula(corr), op=op)
end

function convIndep(x::pbox, y::pbox; op = +)

    if (op == -) return convIndep(x, negate(y), op=+);end
    if (op == /) return convIndep(x, reciprocate(y), op=*);end


    x = makepbox(x);
    y = makepbox(y);

    m = x.n;
    p = y.n;
    n = min(ProbabilityBoundsAnalysis.steps, m*p);
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

    ll = c;	# Probably can remove

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

    bounded = min(x.bounded,y.bounded);

    return pbox(Zu, Zd, ml = ml, mh = mh, vl=vl, vh=vh, bounded = bounded);
    #return ([Zu, Zd, ml, mh, vl, vh, "$(x.dids) $(y.dids)"]);

end


function convPerfect(x::pbox, y::pbox; op = +)

    if (op)∈([-,/])
        cu = map(op, x.u[:],y.d[:]);
        cd = map(op, x.d[:],y.u[:]);
    else
        cu = map(op, x.u[:], y.u[:]);
        cd = map(op, x.d[:], y.d[:]);
    end

    #cu = map(op, x.u[:],y.u[:]);
    #cd = map(op, x.d[:],y.d[:]);
    scu = sort(cu);
    scd = sort(cd);

    # Here we will also need the moment propagation

    bounded = min(x.bounded,y.bounded);

    if (all(cu == scu) && all(cd == scd))
        return pbox(scu, scd, bounded = bounded)
    else return pbox(scu, scd, bounded= bounded)
      end
end

function convOpposite(x::pbox, y::pbox; op = +)


    if (op)∈([-,/])
        cu = map(op, x.u[:],y.d[end:-1:1]);
        cd = map(op, x.d[:],y.u[end:-1:1]);
    else
        cu = map(op, x.u[:], y.u[end:-1:1]);
        cd = map(op, x.d[:], y.d[end:-1:1]);
    end

    #cu = map(op, x.u[:], y.u[end:-1:1]);
    #cd = map(op, x.d[:], y.d[end:-1:1]);
    scu = sort(cu);
    scd = sort(cd);

    # Here we will also need the moment propagation

    bounded = min(x.bounded,y.bounded);

    return pbox(scu, scd, bounded=bounded)
end

function convFrechet(x::pbox, y::pbox; op = +)

    if (op == -) return (convFrechet(x,negate(y),op = +));end
    if (op == /) return (convFrechet(x,reciprocate(y), op = *));end
    if (op == *) if (straddlingzero(x) || straddlingzero(y)) return (imp(balchProd(x,y), convFrechetNaive(x, y, op = *))); end; end
    ## Unsure about the above line. It looks like if it straddles 0, we need to do the naive frechet and the balch prod (?) and impose one on the other

    x = makepbox(x);
    y = makepbox(y);

    zu = zeros(ProbabilityBoundsAnalysis.steps);
    zd = zeros(ProbabilityBoundsAnalysis.steps);

    for i = 1:ProbabilityBoundsAnalysis.steps

        j = i:ProbabilityBoundsAnalysis.steps;
        k = ProbabilityBoundsAnalysis.steps:-1:i;
        zd[i] = minimum(map(op, x.d[j],y.d[k]));

        j = 1:i;
        k = i:-1:1;
        zu[i] = maximum(map(op, x.u[j], y.u[k]));
    end
    #mean

    ml = -Inf;
    mh = Inf;

    if (op)∈([+,-])                 # We should be able to include * /  once we have momement prop
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

    bounded = min(x.bounded,y.bounded);
    return pbox(zu, zd, ml = ml, mh = mh, vl = vl, vh = vh, bounded=bounded);
end



##
#   Can be parred
##
function sigma(x::pbox, y ::pbox; op = +,  C = πCop()::AbstractCopula)

    #if (op == -) return sigma(x, negate(y), op=+, C = C); end
    #if (op == /) return sigma(x, reciprocate(y), op = *, C = C);end

    #x = makepbox(x);
    #y = makepbox(y);

    xus = x.u; xds = x.d;
    yus = y.u; yds = y.d;

    CU = C.cdfU; CD = C.cdfD;

    bounded = min(x.bounded,y.bounded)  # Is the output bounded on each side?

    uCalc = Distributed.@spawn ProbabilityBoundsAnalysis.calcSigmaBoundLeft(xus, yus, op, CU, bounded[1])
    dCalc = Distributed.@spawn ProbabilityBoundsAnalysis.calcSigmaBoundRight(xds, yds, op, CD, bounded[2])

    zu = fetch(uCalc);
    zd = fetch(dCalc);

    # Mean tranforms the same as independence
	ml = -Inf;
    mh = Inf;

    if (op)∈([+, -, *])
        ml = map(op, x.ml, y.ml);
        mh = map(op, x.mh, y.mh);
    end

    # Variance does not
    vl = 0;
    vh = Inf;

    if !(C.cdfU == C.cdfD);
        z1 = deepcopy(zu);  z2 = deepcopy(zd);
        zu[2:end]   = min.(z1[2:end], z2[1:end-1]);
        zd[1:end-1] = max.(z1[2:end], z2[1:end-1]);
    end

    return pbox(zu, zd, ml = ml, mh = mh, vl = vl, vh = vh, bounded = bounded);

end
convCorr(x::pbox, y ::pbox; op = +,  C = πCop()::AbstractCopula) = sigma(x, y, op = op, C = C)

function calcSigmaBoundLeft(x, y , op,  C, bounded)

    xs = x[1:end-1]; ys = y[1:end-1];

    Ns = ProbabilityBoundsAnalysis.steps

    zs = zeros(Ns-1, Ns-1)    # cart prod
    z  = zeros(Ns)            # output vector
    @threads for i = 1:Ns-1
        for j = 1:Ns-1
            zs[i,j] = op(xs[i], ys[j])
        end
    end

    zs = zs[:];

    masses = C[2:end, 2:end] + C[1:end-1, 1:end-1] - C[1:end-1, 2:end] - C[2:end, 1:end-1]

    perm = sortperm(zs)

    zs = zs[perm]
    masses = masses[:][perm]

    z[1] = zs[1]; z[end] = zs[end]

    cdf = cumsum(masses);

    # Shorthand for if, else.
    bounded ? is = ProbabilityBoundsAnalysis.ii() : is = ProbabilityBoundsAnalysis.iii()

    for i = 2:Ns    # Cannot be threaded
        these   = findall(is[i-1] .<= cdf .<= is[i]);
        !isempty(these) ? z[i] = minimum(zs[these]) : z[i] = z[i-1]
    end

    return z

end


function calcSigmaBoundRight(x, y , op,  C, bounded)


    xs = x[2:end]; ys = y[2:end];


    Ns = ProbabilityBoundsAnalysis.steps

    zs = zeros(Ns-1, Ns-1)    # cart prod
    z  = zeros(Ns)            # output vector
    @threads for i = 1:Ns-1
        for j = 1:Ns-1
            zs[i,j] = op(xs[i], ys[j])
        end
    end

    zs = zs[:]
    masses = C[2:end, 2:end] + C[1:end-1, 1:end-1] - C[1:end-1, 2:end] - C[2:end, 1:end-1]

    perm = sortperm(zs)

    zs = zs[perm]
    masses = masses[:][perm]

    z[1] = zs[1]; z[end] = zs[end]

    cdf = reverse(1 .- cumsum(reverse(masses)));

    bounded ? is =  ProbabilityBoundsAnalysis.jj() : is = ProbabilityBoundsAnalysis.jjj()


    for i = 2:Ns    # Cannot be threaded, because of 2nd line
        these   = findall(is[i-1] .<= cdf .<= is[i]);
        !isempty(these) ? z[i-1] = maximum(zs[these]) : (i != 2 ?  z[i-1] = z[i-2] : continue)
    end

    return z

end



function sigmaOld(x::pbox, y ::pbox; op = +,  C = πCop()::AbstractCopula)

    #if (op == -) return sigma(x, negate(y), op=+, C = C); end
    #if (op == /) return sigma(x, reciprocate(y), op = *, C = C);end

    x = makepbox(x);
    y = makepbox(y);

    Ns = ProbabilityBoundsAnalysis.steps

    zus = [map(op, ux, uy) for ux in x.u[1:end-1], uy in y.u[1:end-1]]        # Carteesian products
    zds = [map(op, dx, dy) for dx in x.d[2:end], dy in y.d[2:end]]

    zd = zeros(Ns);
    zu = zeros(Ns);

    zd[1] = minimum(zds); zd[end] = maximum(zds);
    zu[1] = minimum(zus); zu[end] = maximum(zus);

    zus = zus[:];
    zds = zds[:];

    uMasses = C.cdfU[2:end, 2:end] + C.cdfU[1:end-1, 1:end-1] - C.cdfU[1:end-1, 2:end] - C.cdfU[2:end, 1:end-1]
    dMasses = C.cdfD[2:end, 2:end] + C.cdfD[1:end-1, 1:end-1] - C.cdfD[1:end-1, 2:end] - C.cdfD[2:end, 1:end-1]

    pU = sortperm(zus)
    pD = sortperm(zds)

    zus = zus[pU]
    zds = zds[pD]

    uMasses = uMasses[:][pU]
    dMasses = dMasses[:][pD]

    #uCdf = zeros(length(uMasses)+1); dCdf = ones(length(dMasses)+1);

    #push!(zus, zu[end])
    #pushfirst!(zds, zd[1])

    #uCdf[2:end] = cumsum(uMasses);
    #dCdf[1:end-1] = cumsum(dMasses);

    uCdf = cumsum(uMasses);
    dCdf = reverse(1 .- cumsum(reverse(dMasses)));

    #is = range(0, stop = 1, length = Ns); #js = range(0, stop = 1, length = Ns)

    bounded = min(x.bounded,y.bounded);

    is = ProbabilityBoundsAnalysis.iii();
    js = ProbabilityBoundsAnalysis.jjj();
    if bounded[1]; is = ProbabilityBoundsAnalysis.ii(); end
    if bounded[2]; js = ProbabilityBoundsAnalysis.jj(); end

    for i = 2:Ns

        ups   = findall(is[i-1] .<= uCdf .<= is[i]);
        downs = findall(js[i-1] .<= dCdf .<= js[i]);

        if !isempty(ups);   zu[i]   = minimum(zus[ups]);   else zu[i]   = zu[i-1]; end
        if !isempty(downs); zd[i-1] = maximum(zds[downs]); elseif i!=2; zd[i-1] = zd[i-2]; end

    end

    # Mean tranforms the same as independence
	ml = -Inf;
    mh = Inf;

    if (op)∈([+, -, *])
        ml = map(op, x.ml, y.ml);
        mh = map(op, x.mh, y.mh);
    end

    # Variance does not
    vl = 0;
    vh = Inf;

    if !(C.cdfU == C.cdfD);
        z1 = deepcopy(zu);  z2 = deepcopy(zd);
        zu[2:end]   = min.(z1[2:end], z2[1:end-1]);
        zd[1:end-1] = max.(z1[2:end], z2[1:end-1]);
    end

    return pbox(zu, zd, ml = ml, mh = mh, vl = vl, vh = vh, bounded = bounded);

end


function tauRho(x::pbox, y::pbox; op = +, C = W():: AbstractCopula)

    #if (op == -) return (tauRho(x, negate(y), rotate(C), +));end             # Odd behaviour, negating has no effect if we don't do something to copula
    #if (op == /) return (tauRho(x, reciprocate(y), rotate(C), *));end
    #if (op == *) if (straddlingzero(x) || straddlingzero(y)) return (throw(ArgumentError("Not sure if straddles"))); end; end
    ## Unsure about the above line. It looks like if it straddles 0, we need to do the naive frechet and the balch prod (?) and impose one on the other

    x = makepbox(x);
    y = makepbox(y);

    Ns = ProbabilityBoundsAnalysis.steps

    zd = zeros(Ns);
    zu = zeros(Ns);

    zds = [map(op, dx, dy) for dx in x.d, dy in y.d]        # Carteesian products
    zus = [map(op, ux, uy) for ux in x.u, uy in y.u]

    is = ProbabilityBoundsAnalysis.ii(); js = ProbabilityBoundsAnalysis.jj()

    cop = left.(C(js,js))
    dual = right.([is[i] + is[j] - C(is[i], is[j]) for i in 1:Ns, j in 1:Ns])

    cop = cop[:]
    dual = dual[:]

    zds = zds[:]
    zus = zus[:]

    for i = 1:Ns

        downs = findall( cop .>= js[end + 1 - i] );
        ups   = findall( dual .<= is[i] );


        if !isempty(downs);
            zd[i] = minimum(zds[downs]);
            deleteat!(cop, downs);
            deleteat!(zds, downs);
        else zd[i] = zd[i-1]; end


        if !isempty(ups);
            zu[i]   = maximum(zus[ups]);
            deleteat!(dual, ups);
            deleteat!(zus, ups);
        else zu[i]   = zu[i-1]; end



    end

    bounded = min(x.bounded,y.bounded);

    zd = reverse(zd)

    return pbox(zu, zd, bounded = bounded);

end

function tauRhoOld2(x::pbox, y::pbox; op = +, C = W():: AbstractCopula)

    #if (op == -) return (tauRho(x, negate(y), rotate(C), +));end             # Odd behaviour, negating has no effect if we don't do something to copula
    #if (op == /) return (tauRho(x, reciprocate(y), rotate(C), *));end
    #if (op == *) if (straddlingzero(x) || straddlingzero(y)) return (throw(ArgumentError("Not sure if straddles"))); end; end
    ## Unsure about the above line. It looks like if it straddles 0, we need to do the naive frechet and the balch prod (?) and impose one on the other

    x = makepbox(x);
    y = makepbox(y);

    Ns = ProbabilityBoundsAnalysis.steps

    zd = zeros(Ns);
    zu = zeros(Ns);

    zds = [map(op, dx, dy) for dx in x.d, dy in y.d]        # Carteesian products
    zus = [map(op, ux, uy) for ux in x.u, uy in y.u]

    is = range(0, stop = 1, length = Ns); js = range(0, stop = 1, length = Ns)

    cop = C.cdfD;
    dual = [is[i] + js[j] - cop[i,j] for i in 1:Ns, j in 1:Ns]

    downs = findall(cop .== 1); ups = findall(dual .== 0);

    zd[end] = minimum(zds[downs]); zu[1] = maximum(zus[ups])

    for i = 2:Ns

        downs = findall( is[i-1] .< cop .<= is[i]);
        ups   = findall( js[i-1] .<= dual .<= js[i]);

        zd[i-1] = minimum(zds[downs]);
        zu[i]   = maximum(zus[ups]);

    end

    bounded = min(x.bounded,y.bounded);

    return pbox(zu, zd, bounded = bounded);

end


function tauRhoNew(x::pbox, y::pbox; op = +, C = W():: AbstractCopula)

    #if (op == -) return (tauRho(x, negate(y), rotate(C), +));end             # Odd behaviour, negating has no effect if we don't do something to copula
    #if (op == /) return (tauRho(x, reciprocate(y), rotate(C), *));end
    #if (op == *) if (straddlingzero(x) || straddlingzero(y)) return (throw(ArgumentError("Not sure if straddles"))); end; end
    ## Unsure about the above line. It looks like if it straddles 0, we need to do the naive frechet and the balch prod (?) and impose one on the other

    x = makepbox(x);
    y = makepbox(y);

    Ns = ProbabilityBoundsAnalysis.steps

    zd = zeros(Ns);
    zu = zeros(Ns);

    zds = [map(op, dx, dy) for dx in x.d, dy in y.d]        # Carteesian products
    zus = [map(op, ux, uy) for ux in x.u, uy in y.u]

    is = range(0, stop = 1, length = Ns); js = range(0, stop = 1, length = Ns)

    #is = reverse(is)

    cop = C.cdfD;
    dual = [is[i] + js[j] - cop[i,j] for i in 1:Ns, j in 1:Ns]

    cop = cop[:]
    dual = dual[:]

    zds = zds[:]
    zus = zus[:]

    downs = findall(cop .== 1); ups = findall(dual .== 0);

    zd[end] = minimum(zds[downs]); zu[1] = maximum(zus[ups])

    deleteat!(cop, downs)
    deleteat!(dual, ups)

    deleteat!(zds, downs)
    deleteat!(zus, ups)


    for i = 2:Ns

        downs = findall( cop .<= is[i]);
        ups   = findall( dual .< js[i]);

        if !isempty(downs);
            zd[i-1] = minimum(zds[downs]);
            deleteat!(cop, downs);
            deleteat!(zds, downs);
        elseif i!=2; zd[i-1] = zd[i-2]; end

        if !isempty(ups);
            zu[i]   = maximum(zus[ups]);
            deleteat!(dual, ups);
            deleteat!(zus, ups);
        else zu[i]   = zu[i-1]; end

    end

    bounded = min(x.bounded,y.bounded);

    return pbox(zu, zd, bounded = bounded);

end

function tauRho2(x::Real, y::Real, C:: AbstractCopula, op = +)

    #if (op == -) return (tauRho(x,negate(y), C, +));end             # Odd behaviour, negating has no effect if we don't do something to copula
    #if (op == /) return (tauRho(x,reciprocate(y), C, *));end
    #if (op == *) if (straddlingzero(x) || straddlingzero(y)) return (throw(ArgumentError("Not sure if straddles"))); end; end
    ## Unsure about the above line. It looks like if it straddles 0, we need to do the naive frechet and the balch prod (?) and impose one on the other

    x = makepbox(x);
    y = makepbox(y);

    Ns = ProbabilityBoundsAnalysis.steps

    zd = zeros(Ns);
    zu = zeros(Ns);

    zds = [map(op, dx, dy) for dx in x.d, dy in reverse(y.d)]        # Carteesian products
    zus = [map(op, ux, uy) for ux in x.u, uy in reverse(y.u)]

    is = range(0, stop = 1, length = Ns); js = range(0, stop = 1, length = Ns)

    cop = C.cdf;
    dual = [is[i] + js[j] - cop[i,j] for i in 1:Ns, j in 1:Ns]

    downs = findall(cop .== 1); ups = findall(dual .== 0);

    zd[end] = minimum(zds[downs]); zu[1] = maximum(zus[ups])

    for i = 2:Ns

        downs = findall( is[i-1] .<= cop .<= is[i]);
        ups   = findall(js[i-1] .<= dual .<= js[i]);

        zd[i-1] = minimum(zds[downs]);
        zu[i]   = maximum(zus[ups]);

    end

    bounded = min(x.bounded,y.bounded);

    return pbox(zu, zd, bounded=bounded);

end



###
#   Scalars and some univariate functions
###

function shift(x :: pbox, ss :: Real)
     if (x.shape) ∈ (["uniform","normal","cauchy","triangular","skew-normal"]) s = x.shape; else s = ""; end
    return pbox(ss .+ x.u, ss .+ x.d, shape=s, name="", ml = x.ml+ss, mh=x.mh+ss,
    vl=x.vl, vh=x.vh, bounded = x.bounded)
end

function mult(x::pbox, m :: Real)
    if (x.shape) ∈ (["uniform","normal","cauchy","triangular","skew-normal"]) s = x.shape; else s = ""; end
    if ((x.shape) ∈ (["exponential","lognormal"]) && 0 <= x.u[1]) s = x.shape; else s = ""; end
    if (m < 0) return negate(mult(x,abs(m))) end
    return pbox(m*x.u, m*x.d, shape=s, name="", ml=m*x.ml, mh=m*x.mh, vl=(m^2)*x.vl,
    vh=(m^2)*x.vh, bounded = x.bounded)   ################## mean if m<0
end

function negate(x)
    if (ispbox(x))
        if ((x.shape)∈(["uniform", "normal", "cauchy", "triangular"])) s = x.shape; else s = ""; end
        return pbox(-x.d[end:-1:1],-x.u[end:-1:1],shape=s,name = "", ml=-x.mh, mh=-x.ml,
        vl=x.vl, vh=x.vh, bounded = reverse(x.bounded));
    end
    return -x;
end

function complement(x::pbox)
    if ((x.shape)∈(["uniform", "normal", "cauchy", "triangular", "skew-normal"])) s = x.shape; else s = ""; end
    return pbox(1 .-x.d[end:-1:1],1 .-x.u[end:-1:1],shape=s,name = "", ml=1-x.mh,
    mh=1-x.ml, vl=x.vl, vh=x.vh, bounded = x.bounded);
end

function reciprocate(x)
    if ispbox(x)
        if ((x.shape)∈(["Cauchy","{min, max, median}","{min, max, percentile}","{min, max}"]))  sh = x.shape;
        elseif (x.shape == "pareto")    sh = "power function";
        elseif (x.shape == "power function")    sh = "pareto";
        else sh = "";
        end

        #=
        if (left(x) <= 0 && right(x) >= 0)
            return NaN
        else if (left(x)>0)
            myMean = transformMean(x,reciprocate(), false, true);
            myVar = transformVar(x,reciprocate(), false, true);
        else
            myMean = transformMean(x,reciprocate(), false, false);
            myVar = transformVar(x,reciprocate(), false, false);
        end
        =#

        myMean = interval(x.ml, x.mh);
        myVar = interval(x.vl, x.vh);

        return pbox(1 ./reverse(x.d[:]), 1 ./ reverse(x.u[:]), shape = sh, name="", ml=left(myMean),
        mh=right(myMean), vl=left(myVar), vh=right(myVar), bounded = [true, true]);
    end
    return 1/x;
end

function intUnivariate(x :: pbox, op, bounded = [false, false])
    Ints = interval.(x.u,x.d)
    yInts = op.(Ints)

    yu = sort(left.(yInts))
    yd = sort(right.(yInts))

    return pbox(yu, yd, bounded = bounded)
end
#=
function intUnivariate(x :: pbox, op, bounded = [false, false])
    Ints = interval.(x.u[1:end-1],x.d[2:end])
    yInts = op.(Ints)

    xRange = interval(x.u[1], x.d[end])

    yRange = op(xRange)

    yu = sort(left.(yInts))
    yd = sort(right.(yInts))

    return pbox([yRange.lo; yu], [yd; yRange.hi], bounded = bounded)
end
=#
for op in (:sin, :cos, :tan, :sinh, :cosh, :tanh, :asin, :acos, :atan, :exp, :log)
    bounded = [false, false]
    if op ∈ (:sin, :cos); bounded = [true, true]; end
    if op ∈ (:exp, :log); bounded = [true, false]; end
    @eval ($op)(x :: pbox) = intUnivariate(x, $op, $bounded)
end

^(x::pbox, y ::Real) = intUnivariate(x,x->x^y)
^(x::pbox, y ::Integer) = intUnivariate(x,x->x^y)


defaultCorr = interval(-1, 1);

-(x::pbox) = negate(x);
/(x::pbox) = reciprocate(x);

+(x::AbstractPbox, y::AbstractPbox) = conv(x,y, op = +, corr = defaultCorr); # if(x==y) return 2*x; ????
-(x::AbstractPbox, y::AbstractPbox) = conv(x,y, op = -, corr = defaultCorr); # if(x==y) return 0;   ????
*(x::AbstractPbox, y::AbstractPbox) = conv(x,y, op = *, corr = defaultCorr); # if(x==y) return x^2; ????
/(x::AbstractPbox, y::AbstractPbox) = conv(x,y, op =/, corr = defaultCorr); # if(x==y) return 0;   ????

###
#   Conv of pboxes and intervals
###

# Probably will only need shift for + and - with reals
+(x :: AbstractPbox, y :: AbstractInterval) = conv(x, y, op = +, corr = defaultCorr);
+(x :: AbstractInterval, y :: AbstractPbox) = y + x;

-(x :: AbstractPbox, y :: AbstractInterval) = conv(x, y, op =  -, corr = defaultCorr);
-(x :: AbstractInterval, y :: AbstractPbox) = -y + x;

*(x :: AbstractPbox, y :: AbstractInterval) = conv(x, y, op =  *, corr = defaultCorr);
*(x :: AbstractInterval, y :: AbstractPbox) = y * x;

/(x :: AbstractPbox, y :: AbstractInterval) = conv(x, y, op = /, corr = defaultCorr);
/(x :: AbstractInterval, y :: AbstractPbox) = reciprocate(y) * x;

###
#   Conv of pboxes and reals
###

# Probably will only need shift for + and - with reals
+(x :: AbstractPbox, y :: Real) = shift(x, y);
+(x :: Real, y :: AbstractPbox) = y + x;

-(x :: AbstractPbox, y :: Real) = shift(x,-y);
-(x :: Real, y :: AbstractPbox) = -y + x;

*(x :: AbstractPbox, y :: Real) = mult(x,y);
*(x :: Real, y :: AbstractPbox) = y*x;

/(x :: AbstractPbox, y :: Real) = mult(x,1/y);
/(x :: Real, y :: AbstractPbox) = reciprocate(y) * x;



##################################################################################
#
# This (naive) Frechet convolution will work for multiplication of distributions
# that straddle zero.  Note, however, that it is NOT optimal. It can probably be
# improved considerably.  It does use the Frechet moment propagation formulas.
# Unfortunately, moment propagation has not be implemented yet in the R library.
# To get OPTIMAL bounds for the envelope, we probably must use Berleant's linear
# programming solution for this problem.
#
##################################################################################


function convFrechetNaive(x::Real, y::Real; op = *)

    if (op == +) return (convFrechet(x, y, op = +));end
    if (op == -) return (convFrechet(x,negate(y), op = +));end
    if (op == /) return (convFrechetNaive(x,reciprocate(y), op = *));end

    x = makepbox(x);
    y = makepbox(y);
    n = x.n;

    Y = repeat(y.d[:], inner = n);
    X = repeat(x.d[:], outer = n);
    c = sort(map(op,X,Y));
    Zd = c[(n*n - n + 1): n*n];

    Y = repeat(y.u[:], inner = n);
    X = repeat(x.u[:], outer = n);
    c = sort(map(op,X,Y));
    Zu = c[1:n];

    #mean
    m = mean(x) * mean(y);          # Should maybe be op(mean(x),mean(y))
    a = sqrt(var(x) * var(y));      # Simularly
    ml = m - a;
    mh = m + a;

    VK = VKmeanproduct(x,y);
    m = imp(pbox(interval(ml,mh)), VK);

    # Variance
    vl = 0;
    vh = Inf;

    return pbox(Zu, Zd,  ml = left(m), mh = right(m), vl=vl, vh=vh);

end


###
#   Still needed
###

function balchProd(x::pbox, y::pbox)

    if (straddles(x) && straddles(y))
        x0 = left(x);
        y0 = left(y);
        xx0 = x - x0;
        yy0 = y - y0;
        a = convFrechet(xx0, yy0,op= *);
        b = convFrechet(y0 * xx0, x0 * yy0, op = +);
        return convFrechet(a,b, op = +) + x0 * y0;
    end

    if straddles(x)
        x0 = left(x);
        xx0 = x - x0;
        a = convFrechet(xx0, y, op = *);
        b = x0 * y;
        return convFrechet(a, b, op = +);
    end

    if straddles(y)
        y0 = left(y)
        yy0 = y - y0
        a = convFrechet(x, yy0, op = *)
        b = x * y0
        return convFrechet(a, b, op = +)
    end
    return convFrechet(x, y, op = *);
end


function copulaSigma( x :: pbox, y :: pbox; op = +,  C = πCop()::AbstractCopula)

    # Another option is to par the i,j search and return the vector of indicies

    Fz = sigma(x, y, op = op, C = C);

    ux = x.u; uy = y.u; cdf1 = C.cdfU;
    uz = Fz.u;

    dx = x.d; dy = y.d; cdf2 =C.cdfD;
    dz = Fz.d

    proc1 = Distributed.@spawn calcCbound(ux, uy, uz, cdf1, op, true)
    proc2 = Distributed.@spawn calcCbound(dx, dy, dz, cdf2, op, false)

    cdfU = fetch(proc1)
    cdfD = fetch(proc2)

    return copula(cdfU, cdfD), Fz

end

##
# Note in this version of the algorithm everything is a julia base type, ie x == x.u or x.d
##
function calcCbound(x, y, z, C, op, left)

    left ? (xs = x[1:end-1]; ys = y[1:end-1]) : (xs = x[2:end]; ys = y[2:end])

    Nx = length(x); Ny = length(y);
    zus = zeros(Nx-1, Ny-1)
    @threads for i = 1:Nx-1            # Threaded
        for j = 1:Ny-1
            zus[i,j] = op(xs[i], ys[j])
        end
    end

    Ns = size(C)[1]
    cdfC = zeros(Ns, Ns);

    masses = C[2:end, 2:end] + C[1:end-1, 1:end-1] - C[1:end-1, 2:end] - C[2:end, 1:end-1]

    @threads for i = 2:Ns             # Cycle through u's
       for j = 2:Ns                   # Cycle through v's
            indexs = findall((z[i-1] .<= zus .< z[i]));           # This findall may be parallelised
            indexs2 = j-1

            these = getindex.(getproperty.(indexs,:I), 2) .== indexs2       # Quite fast

            us = sum(masses[indexs[these]])

            cdfC[i,j] = us
        end
    end

    return cumsum(cumsum(cdfC, dims=1), dims = 2)
end


function calcCbound2(x, y, z, C, op, left)

    xs = x[2:end]; ys = y[2:end]
    if left; xs = x[1:end-1]; ys = y[1:end-1];end

    #zus =  [map(op, ux, uy) for ux in xs, uy in ys]        # Carteesian products

    Nx = length(x); Ny = length(y);
    zus = zeros(Nx-1, Ny-1)
    for i = 1:Nx-1            # Threaded
        for j = 1:Ny-1
            zus[i,j] = op(xs[i], ys[j])
        end
    end

    Ns = size(C)[1]
    cdfC = zeros(Ns, Ns);

    masses = C[2:end, 2:end] + C[1:end-1, 1:end-1] - C[1:end-1, 2:end] - C[2:end, 1:end-1]

    for i = 2:Ns             # Cycle through u's
       for j = 2:Ns         # Cycle through v's
            indexs = findall((z[i-1] .<= zus .< z[i]));           # This findall may be parallelised
            #indexs2 = findall(x[j-1] .<= x .< x[j]);
            indexs2 = j-1

            these = getindex.(getproperty.(indexs,:I), 2) .== indexs2       # Very fast

            us = sum(masses[indexs[these]])

            cdfC[i,j] = us
            #cdfC[i,j] = us + cdfC[i,j-1] + cdfC[i-1,j] - cdfC[i-1,j-1];
        end
    end

    #cdfC[2:end, 2:end] = [cdfC[i,j] + cdfC[i,j-1] + cdfC[i-1,j] - cdfC[i-1,j-1] for i = 2:Ns, j = 2:Ns]

    return cumsum(cumsum(cdfC, dims=1), dims = 2)
end

## Parallel find
# Not quite working yet
##
function pfind(x1 :: Real, x2 :: Real, M :: Array{Float64,2}, lo :: Integer = 1, hi :: Integer = size(M)[1])

    if hi - lo < size(M)[1]/10      # If it's in 10% of the size, do normal findall
        #these = x1 .<= M[:,hi] .< x2
        #inds = 1:1:size(M)[1]
        #return CartesianIndex.(hi,inds[these[:]])
        return findall(x1 .<= M[:,lo:hi] .<x2)
    end

    mid = (lo+hi)>>>1
    half1 = @spawn pfind(x1, x2, M, lo, mid)
    half2 = @spawn pfind(x1, x2, M, mid+1,hi)

    return [fetch(half1) ; fetch(half2)]

end


function copulaSigmaSlow( x :: pbox, y :: pbox; op = +,  C = πCop()::AbstractCopula)

    Fz = sigma(x, y, op = op, C = C);

    Ns = size(C.cdfD)[1];

    zus = [map(op, ux, uy) for ux in x.u[1:end-1], uy in y.u[1:end-1]]        # Carteesian products
    zds = [map(op, dx, dy) for dx in x.d[2:end], dy in y.d[2:end]]

    cdfU = zeros(Ns, Ns);
    cdfD = zeros(Ns, Ns);

    uMasses = C.cdfU[2:end, 2:end] + C.cdfU[1:end-1, 1:end-1] - C.cdfU[1:end-1, 2:end] - C.cdfU[2:end, 1:end-1]
    dMasses = C.cdfD[2:end, 2:end] + C.cdfD[1:end-1, 1:end-1] - C.cdfD[1:end-1, 2:end] - C.cdfD[2:end, 1:end-1]

    for i = 2:Ns             # Cycle through u's
        for j = 2:Ns         # Cycle through v's
            indexs = findall((Fz.u[i-1] .<= zus .< Fz.u[i]));           # This findall may be parallelised, split domain of integration over threads
            indexs2 = findall(x.u[j-1] .<= x.u .< x.u[j]);              # This is stupid

            YesNo = [indexs[k][1] .== indexs2 for k = 1:length(indexs)];     # Confusing, but checks which elements share the same indexs
            sums = sum.(YesNo);          # The zeros will be removed from the integration

            newIndx = indexs[sums .> 0];

            us = sum(uMasses[newIndx])

            cdfU[i,j] = us + cdfU[i,j-1] + cdfU[i-1,j] - cdfU[i-1,j-1];

            indexs = findall((Fz.d[i-1] .<= zds .< Fz.d[i]));
            indexs2 = findall(x.d[j-1] .<= x.d .< x.d[j]);

            YesNo = [indexs[k][1] .== indexs2 for k = 1:length(indexs)];     # Confusing, but checks which elements share the same indexs
            sums = sum.(YesNo);          # The zeros will be removed from the integration

            newIndx = indexs[sums .> 0];

            ds = sum(dMasses[newIndx])

            cdfD[i,j] = ds + cdfD[i,j-1] + cdfD[i-1,j] - cdfD[i-1,j-1];

        end
        if mod(i,Ns/100 * 20) == 0.0 println("Completed $(i/Ns * 100)%") end
    end

    return copula(cdfU, cdfD), Fz

end


function copulaSigmaHard( x :: pbox, y :: pbox; op = +,  C = πCop()::AbstractCopula)

    Fz = sigma(x, y, op = +, C = C);

    Ns = size(C.cdfD)[1];

    zus = [map(op, ux, uy) for ux in x.u[1:end-1], uy in y.u[1:end-1]]        # Carteesian products
    zds = [map(op, dx, dy) for dx in x.d[2:end], dy in y.d[2:end]]

    cdfU = zeros(Ns, Ns);
    cdfD = zeros(Ns, Ns);

    #zus = zus[:];
    #zds = zds[:];

    uMasses = C.cdfU[2:end, 2:end] + C.cdfU[1:end-1, 1:end-1] - C.cdfU[1:end-1, 2:end] - C.cdfU[2:end, 1:end-1]
    dMasses = C.cdfD[2:end, 2:end] + C.cdfD[1:end-1, 1:end-1] - C.cdfD[1:end-1, 2:end] - C.cdfD[2:end, 1:end-1]

    for i = 1:Ns             # Cycle through u's
        for j = 1:Ns         # Cycle through v's
            indexs = findall((zus .<= Fz.u[i]));
            indexs2 = findall(x.u .<= x.u[j]);

            YesNo = [indexs[k][1] .== indexs2 for k = 1:length(indexs)];     # Confusing, but checks which elements share the same indexs
            sums = sum.(YesNo);          # The zeros will be removed from the integration

            newIndx = indexs[sums .> 0];

            us = sum(uMasses[newIndx])

            cdfU[i,j] = us

        end
        if mod(i,Ns/100 * 20) == 0.0 println("Completed $(i/Ns * 100)%") end
    end

    return copula(cdfU)
end



function copulaUnary(x :: pbox, op)

    Fz = op(x)

    #xInts = interval.(x.u, x.d)
    xRange = range(0, 1,length = Fz.n)
    xInts = cut.(x,xRange)
    zInts = op.(xInts)

    theseFs = cut.(Fz,xRange)
    rightF = right.(theseFs)
    leftF = left.(theseFs)

    CopU = zeros(Fz.n, x.n)
    CopD = zeros(Fz.n, x.n)

    for i = 1:Fz.n
        for j = 1:x.n
            indexU =  left.(zInts[1:j]) .< leftF[i]
            indexD =  right.(zInts[1:j]) .< rightF[i]
            CopU[i, j] = sum(indexU)/Fz.n
            CopD[i, j] = sum(indexD)/Fz.n
        end
    end
    #CopU[1,:].= CopU[:,1] .=  CopD[1,:] .= CopD[:,1] .= 0;      # Set boundary conditions
    #CopU[end,:].= CopU[:,end] .=  CopD[end,:] .= CopD[:,end] .= collect(xRange);
    return copula(CopU,CopD), Fz
end
	
# return the range of a univariate p-box
function Base.range(p::pbox)
    return interval(p.u[1], p.d[end])
end
