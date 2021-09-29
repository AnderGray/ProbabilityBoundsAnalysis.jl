
using IntervalArithmetic

struct pbox{T <: Real}
    u :: Vector{T}    # left discrete inverse
    d :: Vector{T}    # right discrete inverse
    m :: Union{Interval{T}, Missing}  # interval mean
    v :: Union{Interval{T}, Missing}  # interval variance
    shape :: String   # shape
end


using Distributions
const steps = 200

function normal(m :: Interval, std :: Interval)

    ps = range(0, 1, length = steps+1)

    ps_u = ps[1:end-1]
    ps_d = ps[2:end]

    u1 = quantile.(Normal(m.lo, std.lo), ps_u)
    u2 = quantile.(Normal(m.lo, std.hi), ps_u)
    u3 = quantile.(Normal(m.lo, std.hi), ps_u)
    u4 = quantile.(Normal(m.hi, std.hi), ps_u)

    d1 = quantile.(Normal(m.lo, std.lo), ps_d)
    d2 = quantile.(Normal(m.lo, std.hi), ps_d)
    d3 = quantile.(Normal(m.hi, std.lo), ps_d)
    d4 = quantile.(Normal(m.hi, std.hi), ps_d)

    u = min.(u1, u2, u3, u4)
    d = max.(d1, d2, d3, d4)

    return pbox(u, d, m, std^2, "normal")
end

function normal(m :: Interval, std :: Interval)
    ps = range(0, 1, length = steps+1)
    ps_u = ps[1:end-1]; ps_d = ps[2:end]
    u1 = quantile.(Normal(m.lo, std.lo), ps_u);
    d1 = quantile.(Normal(m.hi, std.lo), ps_d)
    u2 = quantile.(Normal(m.lo, std.hi), ps_u);
    d2 = quantile.(Normal(m.hi, std.hi), ps_d)

    u = min.(u1, u2);
    d = max.(d1, d2)
    return pbox(u, d, m, std^2, "normal")
end

function cdf(F :: pbox, x :: Real)

    u = F.u; d = F.d; n = length(u);

    cdf_u = 1 - sum(x .< u)/n;
    cdf_d = 1 - sum(x .<= d)/n;

    return interval(cdf_d, cdf_u)
end

import Base:rand

function cut(F :: pbox, p :: Real)

    u = F.u; d = F.d; n = length(u);
    ps = range(0, 1, length = n+1);

    ind_u = findlast(ps[1:end-1] .<= p)
    ind_d = findfirst(ps[2:end] .>= p)

    return interval(u[ind_u], d[ind_d])
end

# a random sample is a random cut
rand(F :: pbox, N :: Int) = cut.([F], rand(N))

function mass(F :: pbox, x :: Interval)

    cdfHi = cdf(F, x.hi)
    cdfLo = cdf(F, x.lo)  # returns intervals

    ub = cdfHi.hi - cdfLo.lo
    lb = max(0, cdfHi.lo - cdfLo.hi)

    return interval(lb, ub)
end


using MomentArithmetic

function binary(F :: pbox, X :: Interval, op)

  u_z = op.(Interval.(F.u), X)
  d_z = op.(Interval.(F.d), X)

  u_z_ = getfield.(u_z, :lo)  # get right edges
  d_z_ = getfield.(d_z, :hi)  # get left edges

  # for moments
  F_range = interval(F.u[1], F.d[end])
  F_moms = Moments(F.m, F.v, F_range);
  Z_moms = op(F_moms, X);

  return pbox(u_z_, d_z_, Z_moms.mean,
  Z_moms.var, F.shape)
end

using BivariateCopulas
using BivariateCopulas: sample

struct BivPbox{T <: Real}

  X :: pbox{T}  # 1st marginal
  Y :: pbox{T}  # 2nd marginal
  C :: copula   # copula

end

function cdf(F :: BivPbox, x, y)


  cdf_X = cdf(F.X, x)   # returns interval
  cdf_Y = cdf(F.Y, y)

  F_cdf_lo = F.C(cdf_X.lo, cdf_Y.lo)[1]
  F_cdf_hi = F.C(cdf_X.hi, cdf_Y.hi)[1]

  return interval(F_cdf_lo, F_cdf_hi)
end

function rand(F :: BivPbox, N :: Int)

    C_samps = sample(F.C, N)

    X_samps = cut.([F.X], C_samps[:, 1])
    Y_samps = cut.([F.Y], C_samps[:, 2])

    return IntervalBox.(X_samps, Y_samps)

end

function mass(F :: BivPbox, box :: IntervalBox)

    x = box.v[1]; y = box.v[2];  # gets intervals

    F_1 = cdf(F, x.lo, y.lo)
    F_2 = cdf(F, x.hi, y.lo)
    F_3 = cdf(F, x.lo, y.hi)
    F_4 = cdf(F, x.hi, y.hi)

    # using intervals
    mass = F_4 - F_3 - F_2 + F_1;

    return max(0, min(1, mass))
end

function mass(F :: BivPbox, box :: IntervalBox)

x = box.v[1]; y = box.v[2];  # gets intervals

F_1 = cdf(F, x.lo, y.lo)
F_2 = cdf(F, x.hi, y.lo)
F_3 = cdf(F, x.lo, y.hi)
F_4 = cdf(F, x.hi, y.hi)

# using intervals
mass = F_4 - F_3 - F_2 + F_1;

return max(0, min(1, mass))
end

function split(F :: pbox)

    n = length(F.u)

    focal_el = interval.(F.u, F.d)
    masses   = ones(n) ./n

    return focal_el, masses
  end

function makepbox(focal_el :: Vector{Interval{Float64}})

    x_lo = getfield.(focal_el, :lo);
    x_hi = getfield.(focal_el, :hi);

    u = sort(x_lo)
    d = sort(x_hi)

    return pbox(u, d, missing, missing, "")
end


function makepbox( x :: Vector{Interval{T}}, w :: Vector{<:Real}, numel = steps) where T <: Real

    ws = w ./sum(w);

    lefts = getfield.(x, :lo)
    rights = getfield.(x, :hi)

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

    ps = range(0, 1, length = numel+1)
    iis = ps[1:end-1]
    jjs = ps[2:end]

    for p in reverse(iis)
        while true
            if pu[j] <= p
                break
            end
            j = j - 1
        end
        u = [su[j]; u]
    end

    j = 1

    for p in jjs
        while true
            if  p <= pd[j]
                break
            end
            j = j + 1
        end
        d = [d; sd[j]]
    end

    return pbox(u, d, missing, missing, "")


end

function unary(F :: pbox, op)

    focal_el, _ = split(F)
    fe_z = op.(focal_el)

    return makepbox(fe_z)
end

function split(F :: BivPbox)

    X = F.X; Y = F.Y; C = F.C

    nx = length(X.u); ny = length(Y.u);

    px = range(0, 1, length = nx + 1);
    py = range(0, 1, length = ny + 1);

    f_el_X = interval.(X.u, X.d)
    f_el_Y = interval.(Y.u, Y.d)

    masses = C(px[2:end], py[2:end]) - C(px[1:end-1], py[2:end]) - C(px[2:end], py[1:end-1]) + C(px[1:end-1], py[1:end-1])

    f_XY   = [IntervalBox(x,y) for x in f_el_X, y in f_el_Y];

    return f_XY[:], masses[:]
  end


  function sigma(F :: BivPbox, op)

    # Convert pbox to a belief function
    fe, masses = split(F)

    # Perform operation on focal elements
    # with interval arithmetic
    fe_z = [op(f.v[1],f.v[2]) for f in fe]

    # Convert back to p-box
    return makepbox(fe_z, masses)
  end

using PyPlot

function plot(s ::pbox, fill = true; name = missing, col = missing, heading = missing, plotting = true, save = false, alpha = 0.2, fontsize = 12)

    #if (ismissing(name)) name = s.id; end

    col1 = "red"; col2 = "black"; fillcol = "grey"
    if !(ismissing(col)); col1 = col2 = fillcol = col;end

    if !plotting; ioff();end
    if (ismissing(name)); fig = figure(figsize=(10,10)); else fig = figure(name,figsize=(10,10));end

    n = length(s.u)
    ax = fig.add_subplot()
    j = (0:(n-1))/n;

    PyPlot.step([s.u[:];s.u[n];s.d[n]], [j;1;1], color = col1, where = "pre");

    i = (1:(n))/n;
    PyPlot.step([s.u[1];s.d[1];s.d[:]], [0;0;i], color = col2,     where = "post");

    if fill
        Xs, Ylb, Yub = prepFillBounds(s);
        ax.fill_between(Xs, Ylb, Yub, alpha=alpha, color =fillcol)
    end
    if !(ismissing(heading)); title(heading, fontsize = fontsize);end
    xticks(fontsize = fontsize); yticks(fontsize = fontsize)
    xlabel("Distribution range",fontsize = fontsize); ylabel("CDF",fontsize=fontsize);

    if save; savefig("$name.png"); close(fig);end
    ion()
end

plot(s :: Union{<:Real, Interval{<:Real}}, fill = true; name = missing, col = missing, alpha = 0.2, fontsize = 12) = plot(makepbox(s), fill, name = name, col = col, alpha = alpha, fontsize = fontsize)

#plot(s ::pbox, fill = true; name = missing, col = missing, alpha = 0.2, fontsize = 12) = plot(makepbox(s), fill, name=name, col=col, alpha=alpha, fontsize=fontsize)

###
#   Prepares bounds for use with fill between
###
function prepFillBounds(x)

    d = x.d; u = x.u;
    n = length(x.u)

    is = range(0,stop =1 , length = n+1)
    di = is[2:end]; ui = is[1:end-1];

    Xs = sort([d; d; u; u]);
    nums = length(Xs)

    Ylb = zeros(nums,1); Yub = zeros(nums,1);

    for i = 1:2:nums
        #indUb = findfirst(Xs[i] .<= d)
        #indLb = findlast(Xs[i]  .>= u)
        indUb = findlast(Xs[i]  .>= u)
        indLb = findfirst(Xs[i] .<= d)

        if ~isempty(indLb)
            Ylb[i] = ui[indLb]
            if Xs[i] ∈ d
                Ylb[i+1] = di[indLb]
            else
                Ylb[i+1] = ui[indLb]
            end
        else
            Ylb[i] = 1
            Ylb[i+1]=1
        end

        if ~isempty(indUb)
            Yub[i+1] = di[indUb]
            if Xs[i] ∈ u
                Yub[i] = ui[indUb]
            else
                Yub[i] = di[indUb]
            end
        else
            Yub[i] = 0
            Yub[i+1]=0
        end

    end

    return Xs, Ylb[:], Yub[:]
end