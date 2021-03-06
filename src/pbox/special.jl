######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Definition of auxiliary and miscilaneous function
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
#
#   Port of R code pba.r by Scott Ferson and Jason O'Rawe, Applied Biomathematics
#   Origional code available at: https://github.com/ScottFerson/pba.r
######

interval(a::AbstractInterval, b :: AbstractInterval) = hull(a, b)

left(x::Real) = x;
right(x::Real) = x;

left(x::AbstractInterval) = x.lo;
right(x::AbstractInterval) = x.hi;

left(x::AbstractPbox)   = x.u[1];
right(x::AbstractPbox)  = x.d[end];

lefts(x::AbstractPbox)   = x.u
rights(x::AbstractPbox)  = x.d

isinterval(x) = return (typeof(x)<:AbstractInterval ||
(typeof(x) <: AbstractPbox &&  ((x.d[1]==x.d[x.n]) &&  (x.u[1]==x.u[(x.n)]))))

ispbox(x) = return typeof(x)<:AbstractPbox

function isscalar(x)
    if (ispbox(x) && (left(x) == right(x))) return true; end
    if (isinterval(x) && (left(x) == right(x))) return true; end
    if (isa(x, Number) && (length(x) == 1) && !isinterval(x) && !ispbox(x)) return true; end
    return false
end

isvacuous(x::AbstractPbox) = return ((all(x.u .== -Inf)) && (all(x.d .== Inf)));
isvacuous(x::AbstractInterval) = return ((x.lo == -Inf) && (x.hi == Inf));

straddles(x:: Union{AbstractInterval, AbstractPbox}) = return ((left(x)<=0) && (0<=right(x)));   # includes zero
straddlingzero(x:: Union{AbstractInterval, AbstractPbox}) = return ((left(x)<0) && (0<right(x)));   # neglects zero as an endpoint

function degenerate(x :: AbstractInterval)
    if x.lo == x.hi; return x.lo; end
    return x
end

issubset(x :: AbstractVector, y :: IntervalBox{N,T}) where {N,T} = ∈(x,y)
issubset(x :: IntervalBox, y :: Interval) = issubset(x[1],y)
issubset(x :: Interval, y :: IntervalBox) = issubset(x,y[1])

function intersect(x:: Union{Float64,Int64}, y :: AbstractInterval)
    if x ∈ y
        return x
    end
    return ∅
end
intersect(x :: AbstractInterval, y ::Union{Float64,Int64}) = intersect(y,x)


function noNesting(x::Array{<:AbstractInterval})

    N = length(x);
    a=copy(x);
    sort!(a,lt=compareByLo);

    for i =2:N
        if a[i].hi < a[i-1].hi
            return false
        end
    end
    return true
end

function touching(a,b)
    if (isscalar(a) && isscalar(b)) return a == b; end
    if (straddles(a - b)) return true; end
    return false
end


# Sorts item, and sorts follow the same way. Function changes item and follow
function doublequick(item :: Array{<:Real,1}, follow :: Array{<:Real,1}, left = 1, right = length(item))

    i = left;
    j = right;

    x = item[Int(round((left+right)/2))];

    while (i<=j)
        while (item[i] < x && i < right) i+=1;end
        while (x < item[j] && j > left) j -=1;end
        if (i <= j)

			y = item[i];
			item[i] = item[j];
			item[j] = y;

            y = follow[i];
			follow[i] = follow[j];
			follow[j] = y;

			i+=1;
			j-=1;
        end
    end

	if (left < j) doublequick(item, follow, left,j);end
	if (i < right) doublequick(item, follow, i,right);end

end


function Base.show(io::IO, z::pbox)

    digits = 5

    ranglb = round(z.u[1], sigdigits = digits);    rangub = round(z.d[end], sigdigits = digits);
    ml = round(z.ml, sigdigits = digits);   mh = round(z.mh, sigdigits = digits)
    vl = round(z.vl, sigdigits = digits);   vh = round(z.vh, sigdigits = digits)

    if (z.u[1] != z.d[end]) statement1 = "[$ranglb, $rangub]"; else statement1 = "$ranglb";end
    if (z.ml != z.mh) statement2 = "[$ml, $mh]"; else statement2 = "$ml";end
    if (z.vl != z.vh) statement3 = "[$vl, $vh]"; else statement3 = "$vl";end
    print(io, "Pbox: \t $(z.name) ~ $(z.shape) ( range=$statement1, mean=$statement2, var=$statement3)");

end


# Only checks sides and moments, for stricter use ==

#= Not working yet. Requires redefinition of ==

function isequal(x::pbox, y::pbox)
    if (x.u == y.u && x.d == y.d && x.ml == y.ml && x.mh == y.mh && x.vl == y.vl && x.vh == y.vh)
        return true;
    else return false; end
end
isequal(x::pbox, y::AbstractInterval) = return isequal(x,makepbox(y));
isequal(x::AbstractInterval, y::pbox) = return isequal(y,x);

isequal(x::pbox, y::Real) = return isequal(x,makepbox(y));
isequal(x::Real, y::pbox) = return isequal(y,x);

=#
