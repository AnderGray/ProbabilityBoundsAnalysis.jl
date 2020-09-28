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

function conv(x::Real, y::Real, op = +; corr =0) 
    if corr == 0; return convIndep(x,y,op); end
    if corr == 1; return convPerfect(x,y,op);end
    if corr ==-1; return convOpposite(x,y,op);end
    if corr ==interval(-1,1); return convFrechet(x,y,op); end
    
    if isinterval(corr)
        Lower = convCorr(x, y, GauCopula(left(corr)), op);
        Upper = convCorr(x, y, GauCopula(right(corr)), op);
        return env(Lower, Upper)
    end

    return convCorr(x, y, GauCopula(corr), op)
end

function convIndep(x::Real, y::Real, op = +)

    if (op == -) return convIndep(x, negate(y), +);end
    if (op == /) return convIndep(x, reciprocate(y), *);end


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

    return pbox(Zu, Zd, ml = ml, mh = mh, vl=vl, vh=vh, dids="$(x.dids) $(y.dids)",bounded = bounded);
    #return ([Zu, Zd, ml, mh, vl, vh, "$(x.dids) $(y.dids)"]);

end

#=
function sigma(x::Real, y ::Real, C::AbstractCopula, op = +)

    if (op == -) return convCorr(x, negate(y), C, +);end
    if (op == /) return convCorr(x, reciprocate(y), C, *);end

    x = makepbox(x);
    y = makepbox(y);

    Ns = ProbabilityBoundsAnalysis.steps

    zd = zeros(Ns);
    zu = zeros(Ns);

    zds = [map(op, dx, dy) for dx in x.d, dy in y.d]        # Carteesian products
    zus = [map(op, ux, uy) for ux in x.u, uy in y.u]

    uMasses = C.cdfU[2:end,2:end] + C.cdfU[1:end-1,1:end-1] - C.cdfU[1:end-1,2:end] - C.cdfU[2:end,1:end-1]
    dMasses = C.cdfD[2:end,2:end] + C.cdfD[1:end-1,1:end-1] - C.cdfD[1:end-1,2:end] - C.cdfD[2:end,1:end-1]

    is = range(0, stop = 1, length = Ns); js = range(0, stop = 1, length = Ns)

    zd[1] = minimum(zds); zd[end] = maximum(zds);
    zu[1] = minimum(zus); zu[end] = maximum(zus);

    for i = 2:Ns

        downs = findall(is[i-1] .<= cop  .<= is[i]);
        ups   = findall(js[i-1] .<= dual .<= js[i]);
        
        zd[i-1] = minimum(zds[downs]);
        zu[i]   = maximum(zus[ups]);

    end

end
=#


function convCorr(x::Real, y::Real, C:: AbstractCopula, op = +) # This is the same as the conv function except the condensation is different


	if (op == -) return convCorr(x, negate(y), C, +);end
    if (op == /) return convCorr(x, reciprocate(y), C, *);end

    x = makepbox(x);
    y = makepbox(y);

    m = x.n;
    p = y.n;
    n = min(ProbabilityBoundsAnalysis.steps, m*p);
    L = m * p / n;
    c = zeros(m*p);
    Zu = ones(n);
    Zd = ones(n);

	#Probs from copula
    pdf = reshape(C.density,:,1);	# Should be vector of length n*n
	pdf = pdf[:] ./(n*n)					# Gets rid of 2nd array dimension from reshape
	pdfsave = deepcopy(pdf);

    k = 1:n;

    Y = repeat(y.d[:], inner = m);
    for i=1:p
        c[m*(i-1)+1:m*i] = map(op, x.d[:], Y[m*(i-1)+1:m*i]);
    end

	doublequick(c,pdf);
	Zd = condense_d(c, pdf);

	Y = repeat(y.u[:], inner = m);
	c = zeros(length(Y));
	for i=1:p
		c[m*(i-1)+1:m*i] = map(op, x.u[:], Y[m*(i-1)+1:m*i]);
	end

	doublequick(c,pdfsave);
	Zu = condense_u(c, pdfsave);

	# Mean tranforms the same as independence
	ml = -Inf;
    mh = Inf;

    if (op)∈([+,-,*])
        ml = map(op,x.ml,y.ml);
        mh = map(op,x.mh,y.mh);
    end

    # Variance does not
    vl = 0;
    vh = Inf

	# Moment propagation needs to be included here

    bounded = min(x.bounded,y.bounded);

	return pbox(Zu, Zd, ml = ml, mh = mh, vl=vl, vh=vh, dids="$(x.dids) $(y.dids)", bounded = bounded);

end


function condense_d(x :: Array{<:Real,1}, probs :: Array{<:Real,1})

	n = ProbabilityBoundsAnalysis.steps;	d = zeros(n);
	cumulprob = 1.0;	z = n*n;
	laststep = x[z];

	for i = n:-1:1
		pstep = i/n;
		if z > 0
			while ((pstep <= cumulprob) && (z > 0))
				cumulprob = cumulprob - probs[z];			# Can we directly use the cdf?
				z -= 1;
			end
		end
		d[i] = laststep;
		laststep = x[z+1];
	end
	return d
end

function condense_u(x :: Array{<:Real,1}, probs :: Array{<:Real,1})

	n = ProbabilityBoundsAnalysis.steps;	u = zeros(n);
	cumulprob = 0.0; z = 1;
  	laststep = x[z];

	for i = 1:n
		pstep = i/n
		if z <= n*n
			while (pstep >= cumulprob && z<=n*n)
				cumulprob = cumulprob + probs[z];
				z += 1;
			end
		end
		u[i] = laststep;
		laststep = x[z-1];
	end
	return u
end

function convPerfect(x::Real, y::Real, op = +)
    
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
        return pbox(scu, scd,  dids="$(x.dids) $(y.dids) ", bob=x.bob, bounded = bounded)
    else return pbox(scu, scd,  dids="$(x.dids) $(y.dids) ", bounded= bounded)
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
    
    #cu = map(op, x.u[:], y.u[end:-1:1]);
    #cd = map(op, x.d[:], y.d[end:-1:1]);
    scu = sort(cu);
    scd = sort(cd);

    # Here we will also need the moment propagation

    bounded = min(x.bounded,y.bounded);

    return pbox(scu, scd, dids = "$(x.dids) $(y.dids)", bounded=bounded)
end

function convFrechet(x::Real, y::Real, op = +)

    if (op == -) return (convFrechet(x,negate(y),+));end
    if (op == /) return (convFrechet(x,reciprocate(y),*));end
    if (op == *) if (straddlingzero(x) || straddlingzero(y)) return (imp(balchProd(x,y), convFrechetNaive(x, y, *))); end; end
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
    return pbox(zu, zd, ml = ml, mh = mh, vl = vl, vh = vh, dids = "$(x.dids) $(y.dids) ", bounded=bounded);
end


function tauRho(x::Real, y::Real, C:: AbstractCopula, op = +)

    #if (op == -) return (tauRho(x,negate(y), C, +));end             # Odd behaviour, negating has no effect if we don't do something to copula
    #if (op == /) return (tauRho(x,reciprocate(y), C, *));end
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
    vl=x.vl, vh=x.vh, dids=x.dids, bob=perfectdep(x),bounded = x.bounded)
end

function mult(x::pbox, m :: Real)
    if (x.shape) ∈ (["uniform","normal","cauchy","triangular","skew-normal"]) s = x.shape; else s = ""; end
    if ((x.shape) ∈ (["exponential","lognormal"]) && 0 <= x.u[1]) s = x.shape; else s = ""; end
    if (m < 0) return negate(mult(x,abs(m))) end
    return pbox(m*x.u, m*x.d, shape=s, name="", ml=m*x.ml, mh=m*x.mh, vl=(m^2)*x.vl, 
    vh=(m^2)*x.vh, dids=x.dids, bob=perfectopposite(m,x), bounded = x.bounded)   ################## mean if m<0
end

function negate(x)
    if (ispbox(x))
        if ((x.shape)∈(["uniform", "normal", "cauchy", "triangular"])) s = x.shape; else s = ""; end
        return pbox(-x.d[end:-1:1],-x.u[end:-1:1],shape=s,name = "", ml=-x.mh, mh=-x.ml, 
        vl=x.vl, vh=x.vh, dids=x.dids, bob=oppositedep(x), bounded = reverse(x.bounded));
    end
    return -x;
end

function complement(x::pbox)
    if ((x.shape)∈(["uniform", "normal", "cauchy", "triangular", "skew-normal"])) s = x.shape; else s = ""; end
    return pbox(1 .-x.d[end:-1:1],1 .-x.u[end:-1:1],shape=s,name = "", ml=1-x.mh, 
    mh=1-x.ml, vl=x.vl, vh=x.vh, dids=x.dids, bob=oppositedep(x), bounded = x.bounded);
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
        mh=right(myMean), vl=left(myVar), vh=right(myVar), dids=x.dids, bob=oppositedep(x), bounded = [true, true]);
    end
    return 1/x;
end


defaultCorr = 0;

-(x::pbox) = negate(x);
/(x::pbox) = reciprocate(x);

+(x::AbstractPbox, y::AbstractPbox) = conv(x,y,+, corr = defaultCorr); # if(x==y) return 2*x; ????
-(x::AbstractPbox, y::AbstractPbox) = conv(x,y,-, corr = defaultCorr); # if(x==y) return 0;   ????
*(x::AbstractPbox, y::AbstractPbox) = conv(x,y,*, corr = defaultCorr); # if(x==y) return x^2; ????
/(x::AbstractPbox, y::AbstractPbox) = conv(x,y,/, corr = defaultCorr); # if(x==y) return 0;   ????

###
#   Conv of pboxes and intervals
###

# Probably will only need shift for + and - with reals
+(x :: AbstractPbox, y :: AbstractInterval) = conv(x, y, +, corr = defaultCorr);
+(x :: AbstractInterval, y :: AbstractPbox) = y + x;

-(x :: AbstractPbox, y :: AbstractInterval) = conv(x, y, -, corr = defaultCorr);
-(x :: AbstractInterval, y :: AbstractPbox) = -y + x;

*(x :: AbstractPbox, y :: AbstractInterval) = conv(x, y, *, corr = defaultCorr);
*(x :: AbstractInterval, y :: AbstractPbox) = y * x;

/(x :: AbstractPbox, y :: AbstractInterval) = conv(x, y, /, corr = defaultCorr);
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


oppositedep(x::pbox) = -x.bob;
perfectdep(x::pbox) = x.bob;
perfectopposite(m, x::pbox) = if (m<0) return oppositedep(x); else return perfectdep(x);end


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


function convFrechetNaive(x::Real, y::Real, op = *)

    if (op == +) return (convFrechet(x, y,+));end
    if (op == -) return (convFrechet(x,negate(y),+));end
    if (op == /) return (convFrechetNaive(x,reciprocate(y),*));end

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

    return pbox(Zu, Zd,  ml = left(m), mh = right(m), vl=vl, vh=vh, dids="$(x.dids) $(y.dids)");

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
        a = convFrechet(xx0, yy0, *);
        b = convFrechet(y0 * xx0, x0 * yy0, +);
        return convFrechet(a,b,+) + x0 * y0;
    end

    if straddles(x)
        x0 = left(x);
        xx0 = x - x0;
        a = convFrechet(xx0, y, *);
        b = x0 * y;
        return convFrechet(a,b,+);
    end

    if straddles(y)
        y0 = left(y)
        yy0 = y - y0
        a = convFrechet(x, yy0, *)
        b = x * y0
        return convFrechet(a,b,+)
    end
    return convFrechet(x,y,*);
end
