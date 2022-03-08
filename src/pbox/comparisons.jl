######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Adds comparision functions between p-boxes
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
#
######

###
#   Scalar and p-box comparisons
###

function <(x::pbox, y::Real)
    if x.d[end] < y; return true; end
    if y < x.u[1]; return false; end
    return cdf(x, y)
end

function >(x::pbox, y::Real)
    if y < x.u[1]; return true; end
    if x.d[end] < y; return false; end
    return 1 - cdf(x, y)
end

function <=(x::pbox, y::Real)
    if x.d[end] <= y; return true; end
    if y <= x.u[1]; return false; end
    return cdf(x, y)
end

function >=(x::pbox, y::Real)
    if y <= x.u[1]; return true; end
    if x.d[end] <= y; return false; end
    return 1 - cdf(x, y)
end

<(x::Real, y::pbox) = y > x
>(x::Real, y::pbox) = y < x
<=(x::Real, y::pbox) = y >= x
>=(x::Real, y::pbox) = y <= x

###
#   Interval and p-box comparisons
###

function <(x::pbox, y::Interval{T}) where T
    if x.d[end] < y.lo; return true; end
    if y.hi < x.u[1]; return false; end
    z = x - y
    return cdf(z, 0)
end

function <=(x::pbox, y::Interval{T}) where T
    if x.d[end] <= y.lo; return true; end
    if y.hi < x.u[1]; return false; end
    z = x - y
    return cdf(z, 0)
end

function >=(x::pbox, y::Interval{T}) where T
    if y.hi <= x.u[1]; return true; end
    if x.d[end] <= y.lo; return false; end
    z = y - x
    return cdf(z, 0)
end

###
#   p-box and p-box comparisons
###

function <(x::pbox, y::pbox; corr = interval(-1,1))

    if x.d[end] < y.u[1]; return true; end
    if y.d[end] < x.u[1]; return false; end

    z = conv(x, y, op = -, corr = corr)

    val = cdf(z, 0);

    if val == 1; return true; end
    if val == 1..1; return true; end
    if val == 0; return false; end
    if val == 0..0; return false; end

    return val
end

function <=(x::pbox, y::pbox; corr = interval(-1,1))

    if x.d[end] <= y.u[1]; return true; end
    if y.d[end] <= x.u[1]; return false; end

    z = conv(x, y, op = -, corr = corr)

    val = cdf(z, 0);

    if val == 1; return true; end
    if val == 1..1; return true; end
    if val == 0; return false; end
    if val == 0..0; return false; end

    return val
end

>(x::pbox, y::pbox; corr = interval(-1,1)) = <(y, x, corr = corr)
>=(x::pbox, y::pbox; corr = interval(-1,1)) = <(y, x, corr = corr)
