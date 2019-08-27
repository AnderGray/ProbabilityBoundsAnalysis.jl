######
# This file is part of the pba.jl package.
#
# Definition of auxiliary and miscilaneous function
#
#           University of Liverpool
######

left(x::Real) = x;
right(x::Real) = x;

left(x::AbstractInterval) = x.lo;
right(x::AbstractInterval) = x.hi;

left(x::AbstractPbox)   = x.u[1];
right(x::AbstractPbox)  = x.d[end];


isinterval(x) = return (typeof(x)<:AbstractInterval ||
(typeof(x) <: AbstractPbox &&  ((x.d[1]==x.d[x.n]) &&  (x.u[1]==x.u[(x.n)]))))

ispbox(x) = return typeof(x)<:AbstractPbox

isvacuous(x::AbstractPbox) = return ((all(x.u .== -Inf)) && (all(x.d .== Inf)));
isvacuous(x::AbstractInterval) = return ((x.lo == -Inf) && (x.hi == Inf));

straddles(x:: Union{AbstractInterval, AbstractPbox}) = return ((left(x)<=0) && (0>=right(x)));   # includes zero
straddlingzero(x:: Union{AbstractInterval, AbstractPbox}) = return ((left(x)<0) && (0>right(x)));   # neglects zero as an endpoint




function Base.show(io::IO, z::pbox)

    if (z.u[1] != z.d[end]) statement1 = "[$(z.u[1]),$(z.d[end])]"; else statement1 = "$(z.u[1])";end
    if (z.ml != z.mh) statement2 = "[$(z.ml),$(z.mh)]"; else statement2 = "$(z.ml)";end
    if (z.vl != z.vh) statement3 = "[$(z.vl),$(z.vh)]"; else statement3 = "$(z.vl)";end
    print(io, "Pbox: \t $(z.name) ~ $(z.shape) ( range=$statement1, mean=$statement2, var=$statement3)");

end


# Only checks sides and moments, for stricter use ==

#= Not working yet. Requires definition of ==

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
