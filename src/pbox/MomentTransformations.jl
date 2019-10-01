######
# This file is part of the pba.jl package.
#
# Definition of arithmetic and transformations of statistical moments
#
#           University of Liverpool
######


###
# Todo
###


function transformMean(x :: AbstractPbox, t :: Function, posDir :: Bool, pos2Dir :: Bool)
        #Transformation  of the Mean
end

function transformVar(x :: AbstractPbox, t :: Function, posDir :: Bool, pos2Dir :: Bool)
        #Transformation  of the Variance
end



# THE FOLLOWING VLADIK MOMENT ROUTINES SUPPORT NAIVEFRECHETCONVPBOX

#///////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////
#// Vladik's mean formulations
#//
#// Ferson, S. L. Ginzburg, V. Kreinovich and J. Lopex. 2002. Absolute bounds on the mean of sum,
#// product, max, and min: a probabilistic extension of interval arithmetic.
#//
#// Vladik:  On page 19, I don't understand the nature of the two kinds of conditions:  1 contained in
#// the sum p1+p2 (which I call pa+pb), and p1 (pa) overlapping p2 (pb).  Are the implementations ok?
#// Also, I believe that the minuends in the denominators of p_i at the top of page 19 should be
#// right(x_i) rather than right(E_i). On page 21, what is the function $\underbar{f}^u_{min}(p_1,p_2)$ ?
#///////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////

function VKmeanlo(p :: Float64, q :: Float64, k1 :: Float64, k2 :: Float64, k3 :: Float64, k4 :: Float64)
        return max(p + q - 1, 0.0) * k1 + min(p, 1 - q) * k2 + min(1 - p, q) * k3 + max(1 - p - q, 0.0) * k4;
end

function VKmeanup(p :: Float64, q :: Float64, k1 :: Float64, k2 :: Float64, k3 :: Float64, k4 :: Float64)
        return min(p, q) * k1 + max(p - q, 0.0) * k2 + max(q - p, 0.0) * k3 + min(1 - p, 1 - q) * k4;
end

function VKmeanlower(pa :: AbstractInterval, pb :: AbstractInterval, k1 :: Float64, k2 :: Float64, k3 :: Float64, k4 :: Float64)

        lpa = left(pa);
        rpa = right(pa);
        lpb = left(pb);
        rpb = right(pb);
        touch = touching(pa + pb, 1.0);
        p1 = lpa;               p2 = lpb;               ec1 =          VKmeanlo(p1,p2,k1,k2,k3,k4);
        p1 = lpa;               p2 = rpb;               ec1 = min(ec1, VKmeanlo(p1,p2,k1,k2,k3,k4));
        p1 = rpa;               p2 = lpb;               ec1 = min(ec1, VKmeanlo(p1,p2,k1,k2,k3,k4));
        p1 = rpa;               p2 = rpb;               ec1 = min(ec1, VKmeanlo(p1,p2,k1,k2,k3,k4));
        p1 = max(lpa, 1 - rpb); p2 = 1 - p1; if (touch) ec1 = min(ec1, VKmeanlo(p1,p2,k1,k2,k3,k4)); end
        p1 = min(rpa, 1 - lpb); p2 = 1 - p1; if (touch) ec1 = min(ec1, VKmeanlo(p1,p2,k1,k2,k3,k4)); end
        return ec1
end

function VKmeanupper(pa :: AbstractInterval, pb :: AbstractInterval, k1 :: Float64, k2 :: Float64, k3 :: Float64, k4 :: Float64)

        lpa = left(pa);
        rpa = right(pa);
        lpb = left(pb);
        rpb = right(pb);
        touch = touching(pa, pb);
        p1 = lpa;               p2 = lpb;              ec2 =          VKmeanup(p1,p2,k1,k2,k3,k4);
        p1 = lpa;               p2 = rpb;              ec2 = max(ec2, VKmeanup(p1,p2,k1,k2,k3,k4));
        p1 = rpa;               p2 = lpb;              ec2 = max(ec2, VKmeanup(p1,p2,k1,k2,k3,k4));
        p1 = rpa;               p2 = rpb;              ec2 = max(ec2, VKmeanup(p1,p2,k1,k2,k3,k4));
        p1 = max(lpa, 1 - rpb); p2 = 1-p1;  if (touch) ec2 = max(ec2, VKmeanup(p1,p2,k1,k2,k3,k4)); end
        p1 = min(rpa, 1 - lpb); p2 = 1-p1;  if (touch) ec2 = max(ec2, VKmeanup(p1,p2,k1,k2,k3,k4)); end
        return ec2
end

function VKmeanproduct(a, b)

        # Interval ea,eb,ec,pa,pb;
        # double la,ra,lb,rb,k1,k2,k3,k4,ec1,ec2;

        ec = interval(-Inf,Inf);
        ea = mean(a);
        eb = mean(b);
        la = left(a);
        ra = right(a);
        lb = left(b);
        rb = right(b);
        k1 = ra * rb;
        k2 = ra * lb;
        k3 = la * rb;
        k4 = la * lb;
        pa = interval( (left(ea) - la) / (ra - la), (right(ea) - la) / (ra - la) );
        pb = interval( (left(eb) - lb) / (rb - lb), (right(eb) - lb) / (rb - lb) );
        ec = env(VKmeanlower(pa,pb,k1,k2,k3,k4), VKmeanupper(pa,pb,k1,k2,k3,k4));
        return ec
end
