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


function no_nesting(x::Array{<:AbstractInterval})

    N = length(x);
    a=copy(x);
    sort!(a,lt= (x,y) -> (x.lo <= y.lo));

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
