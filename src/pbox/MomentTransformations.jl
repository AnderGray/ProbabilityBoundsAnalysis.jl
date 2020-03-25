######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#	Definition of arithmetic and transformations of statistical moments
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
#
#   Port of R code pba.r by Scott Ferson and Jason O'Rawe, Applied Biomathematics
#   Origional code available at: https://github.com/ScottFerson/pba.r
######

###
#	->	If not included: When the mean of an interval is also included, the bounds on the variance change as:
#		Lower = 0;
#		upper = (µ - lo)(µ - hi)
#
###

##
#	I think I need to look over row's paper to understand this
##
function transformMean(x :: AbstractPbox, t :: Function, posDir :: Bool, pos2Dir :: Bool)
        #Transformation  of the Mean
        if (ismissing(pos2Dir))
                return transformedMeanMissing2ndDeriv(x, t, posDir);
        elseif (pos2Dir)
                return transformedMeanPos2ndDeriv(x, t, posDir);
        else
                newT(a) = (-1) * t(a);
                tempResult = transformedMeanPos2ndDeriv(x, newT, !posderiv)

                if (!ismissing(tempResult))
                        transformedMean_lower = (-1) * (tempResult.hi)
                        transformedMean_upper = (-1) * (tempResult.lo)
                        return interval(transformedMean_lower, transformedMean_upper)
                else
                        return missing
                end
        end
end

function transformedMeanMissing2ndDeriv(x :: AbstractPbox, t :: Function, posDir :: Bool)

        transformedMean_lower = missing;
        transformedMean_upper = missing;

        transformedSum = 0;
        for i in 1:x.n
		transformedSum = transformedSum + t(x.u[i])
	end

	if (posDir)
		transformedMean_lower = transformedSum/(x.n)
	else
		transformedMean_upper = transformedSum/(x.n)
	end

	transformedSum = 0;
        for i in 1:x.n
		transformedSum = transformedSum + t(x.d[i])
	end

	if (posDir)
		transformedMean_upper = transformedSum/(x.n)
	else
		transformedMean_lower = transformedSum/(x.n)
	end

	return interval(transformedMean_lower, transformedMean_upper);
end

function transformedMeanPos2ndDeriv(x :: AbstractPbox, t :: Function, posDir :: Bool)

	transformedMean_lower = missing;
	transformedMean_upper = missing;
	mu_lower = x.ml;
	mu_upper = x.mh;

	sum = 0;
	transformedSum = 0;
	for i in 1:x.n
		sum = sum + x.u[i];
		transformedSum = transformedSum + t(x.u[i])
	end

	computedMu_lower = sum/x.n;
	if (posderiv)
		transformedMean_lower = transformedSum/(x.n)
	else
		transformedMean_upper = transformedSum/(x.n)
	end

	sum = 0;
	transformedSum = 0;
	for i in 1:x.n
		sum = sum + x.d[i];
		transformedSum = transformedSum + t(x.d[i])
	end

	computedMu_upper = sum/x.n;
	if (posderiv)
		transformedMean_upper = transformedSum/(x.n)
	else
		transformedMean_lower = transformedSum/(x.n)
	end

	if (computedMu_lower > mu_upper || computedMu_upper < mu_lower)
		return missing
	elseif (computedMu_lower >= mu_lower && computedMu_upper <= mu_upper)
		return interval(transformedMean_lower,transformedMean_upper)
	else
		if posDir
			if computedMu_lower < mu_lower
				transformedMean_lower = missing;
			end
			if computedMu_upper > mu_upper
				transformedMean_upper = missing;
			end
		else
			if computedMu_lower < mu_lower
				transformedMean_upper = missing;
			end
			if computedMu_upper > mu_upper
				transformedMean_lower = missing;
			end
		end
	end
#=
	if ismissing(transformedMean_lower)

		if posDir mu = mu_lower; else mu = mu_upper; end

		endingPoints = getBoundEndingPoints(x);

		sum = 0;
		transformedSum = 0;

		for i in 1:x.n
			sum = sum + x.u[i];
			transformedSum = transformedSum + t(x.u[i])
		end

		num_middle = 0

		computedMu = sum/(x.n);

		if ( computedMu >= mu_lower && computedMu <= mu_upper)
			transformedMean_lower = transformedSum/(x.n)
		end

		for i in 1:(2*(B.n))

			if (endingPoints[2, i] < 0)
				sum = sum - endingPoints[1, i]
				transformedSum = transformedSum - t(endingPoints[1, i])
				num_middle = num_middle + 1
			else
				sum = sum + endingPoints[1, i]
				transformedSum = transformedSum + t(endingPoints[1, i])
				num_middle = num_middle - 1
			end


			if (num_middle == 0)

				computedMu = sum/(x.n)
				if ( computedMu >= mu_lower && computedMu <= mu_upper)
					currentTransformedMean = transformedSum/(x.n)

					if (ismissing(transformedMean_lower))
						transformedMean_lower = currentTransformedMean
					else
						transformedMean_lower = min(transformedMean_lower, currentTransformedMean)
					end
				else
					value_middle = (mu * (x.n) - sum ) / num_middle

					if ( value_middle >= endingPoints[1, i] &&  value_middle <= endingPoints[1, i+1])
							currentTransformedMean = (transformedSum + t(value_middle)*num_middle)/(x.n)

							if (ismissing(transformedMean_lower))
								transformedMean_lower = currentTransformedMean
							else
								transformedMean_lower = min(transformedMean_lower, currentTransformedMean)
							end
						end
					end
				end
			end
		end
	end

	if ismissing(transformedMean_upper)

		if posDir
			mu = mu_upper;
		else
			mu = mu_lower;
		end
		sum = mu * x.n;

		sum_lower = 0;
		transformedSum_lower = 0;

		for i in 1:(x.n)
			sum_lower = sum_lower + x.d[i]
			transformedSum_lower = transformedSum_lower + t(x.d[i])
    	end

		for i in 1:(x.n)
  			sum_upper = sum_lower
			transformedSum_upper = transformedSum_lower
			sum_lower = sum_lower +  x.u[i] - x.d[i]
			transformedSum_lower = transformedSum_lower + t(x.u[i]) - t(x.d[i])

			if ( sum >= sum_lower &&  sum <= sum_upper)

				a = (sum_upper - sum) / (x.d[i] - x.u[i])

				currentTransformedMean = (transformedSum_upper - a * ( t(x.d[i]) - t(x.u[i]) ) )/(x.n)

				if (ismissing(transformedMean_upper))
					transformedMean_upper = currentTransformedMean
				else
					transformedMean_upper = max(transformedMean_upper, currentTransformedMean)
				end
			end
		end
	end
	=#
	return interval(transformedMean_lower,transformedMean_upper);
end

function transformVar(x :: AbstractPbox, t :: Function, posDir :: Bool, transformedMean :: AbstractInterval)

	if ismissing(pos2Dir)
		return transformedVarianceByIntervalComp(x, t, posDir, transformedMean);
	elseif pos2Dir
		return transformedVariancePos2ndDeriv(x, t, posDir, transformedMean);
	else
		newT(a) = (-1) * t(a);
		newTransformedMean = interval((-1) * (transformedMean.hi), (-1) * (transformedMean.lo));
		return transformedVariancePos2ndDeriv(x, newT, !posDir, newTransformedMean)
	end
end

function transformedVariancePos2ndDeriv(x::AbstractPbox, t :: Function, posDir :: Bool, transformedMean :: AbstractInterval)
end

function transformedVarianceByIntervalComp(x::AbstractPbox, t :: Function, posDir :: Bool, transformedMean :: AbstractInterval)
end

function transformedVarianceByRowesPos2ndDeriv(x::AbstractPbox, t :: Function, posDir :: Bool, transformedMean :: AbstractInterval)
end

function getBoundEndingPoints(x:: AbstractPbox)

	endingPoints = [];
	i_u = 1;
	i_d = 1;

	for i=1:(2*x.n)
		if (i_u > x.n)

			push!(endingPoints,x.d[i_d]);
			push!(endingPoints,1);
			i_d = i_d + 1;

		elseif (x.u[i_u] < x.d[i_d])

			push!(endingPoints,x.u[i_u]);
			push!(endingPoints,-1);
			i_u = i_u + 1;

		else

			push!(endingPoints,x.d[i_d]);
			push!(endingPoints,1);
			i_d = i_d + 1;

		end
	end

	endingPoints = reshape(endingPoints, 2, (2*x.n));
	return endingPoints;

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
