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

function <(x::pbox, y::pbox; corr = interval(-1,1))
    z = conv(x, y, op = -, corr = corr)
    return cdf(z, 0)
end

<=(x::pbox, y::pbox; corr = interval(-1,1)) = <(x::pbox, y::pbox, corr = corr)
>(x::pbox, y::pbox; corr = corr) = <(y, x, corr = corr)
>=(x::pbox, y::pbox; corr = corr) = <(y, x, corr = corr)


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
