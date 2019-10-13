######
# This file is part of the pba.jl package
#
#   Definition of pbox class with constructor methdos
#
#           University of Liverpool
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
#
#   Port of C++ by Scott Ferson (2007)
#   Origional code available at:
#   Based on the paper: "Experimental Uncertainty Estimation and Statistics for
#                               Data Having Interval Uncertainty"
#
######

using IntervalArithmetic, Statistics


# For saving upper or lower bound information
struct BoundItem
    value :: Real
    bound :: Int64      # 0 -> lowerbound, 1 -> upperbound

    function BoundItem( _value :: Real, _bound ::Int64)
        new(_value,_bound)
    end
end

compareByValue(a :: BoundItem, b :: BoundItem) = (return a.value < b.value);
compareByLo(a::AbstractInterval, b::AbstractInterval) = (return a.lo < b.lo);
#isless(a::BoundItem, b::BoundItem) = isless(a.value, b.value) #May be faster to overload the less than

Variance(a::Array{<:AbstractInterval}) = interval(LowerBoundVariance(a),UpperBoundVariance(a));
var(a::Array{<:AbstractInterval}) = Variance(a);

StandardDeviation(a::Array{<:AbstractInterval}) = interval(sqrt(LowerBoundVariance(a)),sqrt(UpperBoundVariance(a)));
std(a::Array{<:AbstractInterval}) = StandardDeviation(a)


function SampleVariance(a::Array{<:AbstractInterval})
    N = length(a);
    q = N/(N-1.0);
    v = Variance(a)
    return interval(q*v.lo, q*v.hi);
end

SampleVar(a::Array{<:AbstractInterval}) = SampleVariance(a);

function SampleStandardDeviation(a::Array{<:AbstractInterval})
    N = length(a);
    q = N/(N-1.0);
    s = StandardDeviation(a);

    return interval(q*s.lo,q*s.hi)
end

SampleStd(a::Array{<:AbstractInterval}) = SampleStandardDeviation(a);

function mean(x::Array{<:AbstractInterval})
    N = length(x);
    sumHi = 0.0; sumLo = 0.0;
    for i = 1:N
        sumLo += x[i].lo
        sumHi += x[i].hi
    end
    return interval(sumLo/N,sumHi/N)
end

function LowerBoundVariance(x :: Array{<:AbstractInterval})

    N = length(x);
    Xm = mean(x);
    Y = [ BoundItem(x[1].lo,0); BoundItem(x[1].hi,1) ];

    for i = 2:N
        push!(Y,BoundItem(x[i].lo,0));
        push!(Y,BoundItem(x[i].hi,1));
    end

    sort!(Y, lt = compareByValue);
    first = true; Nk = 0; Sk = 0.0; Mk = 0.0;
    V = Inf;

    for k = 1:length(Y)-1
        Yk_k1 = interval(Y[k].value,Y[k+1].value);
        if ( !isdisjoint(Xm, Yk_k1) )
            if (first)
                # all upper bound Yi before and include Yk
                for i=1:k
                    #  if upper bound
                    if (Y[i].bound == 1)
                        v = Y[i].value;
                        Nk += 1;
                        Sk += v;
                        Mk += v*v/N;
                    end
                end
                # all lower bound Yi after and include Yk
                for j= (k+1):length(Y)-1    # Last value cannot be a lowerbound
                    # if lower bound
                    if (Y[j].bound == 0)
                        v = Y[j].value;
                        Nk +=1;
                        Sk += v;
                        Mk += v*v/N;
                    end
                end
            else
                # if lower bound, gain a new value (Lose?)
                if (Y[k].bound == 0)
                    v = Y[k].value;
                    Nk -= 1;
                    Sk -= v;
                    Mk -= v*v/N;
                end
                # if upper bound, lose a new value (Gain?)
                if (Y[k].bound == 1)
                    v = Y[k].value;
                    Nk += 1;
                    Sk += v;
                    Mk += v*v/N;
                end
            end
            if (Nk == 0) return 0; end
            ratio = Sk / Nk;
            if (in(ratio, Yk_k1) && in(ratio, Xm))
                Vk = Mk - (Sk * Sk)/(N * Nk);
                if (Vk < V) V = Vk; end
            end
            first = false;
        end
    end
    return V
end

function UpperBoundVariance( x :: Array{<:AbstractInterval})

    N = length(x);
    sort!(x, lt = compareByLo);

    # Check for no nesting, otherwise return -1 (could also return general function)
    if (!noNesting(x)) return SlowUpperBoundVariance(x) end
    M = 0.0; E = 0.0;
    for i = 1:N
        E += x[i].hi / N;
        M += x[i].hi * x[i].hi / N;
    end
    V = M - E * E;
    for k = 1:N
        iN = x[k].lo;
        oUt = x[k].hi;
        M += iN * iN / 3.0 - oUt * oUt / 3.0;
        E += iN / 3.0 - oUt / 3.0;
        v = M - E * E;
        if V < v
            V = v;
        end
    end
    return V;
end

function SlowUpperBoundVariance(a::Array{<:AbstractInterval})

    bounds = setProduct(a)
    V = 0
    for i = 1:length(bounds)
        v = var(bounds[i])
        if (V<v)
            V=v;
        end
    end
    return V;
end



function setProduct(a::Array{<:AbstractInterval})

    N = length(a);
    cardinality = 2^N;
    temp = zeros(cardinality,N);
    for i = 1:N
        temp[:,i] = repeat([a[i].lo,a[i].hi],inner = 2^(i-1), outer = 2^(N-i));
    end
    result = [temp[i,:] for i = 1:cardinality];
    return result;
end


##
#   This is the general algorithm, which is very very slow
#   My own impementation from Scotts sandia report.
#
#   The bellow code is most likely to be removed now
##

function IntervalStatistics(samples)

    if isequal(samples[:,1],samples[:,2])
        means = mean(samples[:,1]);
        vars = var(samples[:,1]);
        intervalMean = [means, means];
        intervalVariance = [vars,vars];
        return [intervalMean, intervalVariance];
     end

     intervalMean = [mean(samples[:,1]),mean(samples[:,2])];
     intervalVariance = [0.0,0.0];

     samplesSorted = zeros(size(samples));

     samplesSorted[:,1] = sort(samples[:,1]);
     samplesSorted[:,2] = sort(samples[:,2]);

     Ys = reshape(samplesSorted',length(samples),1);

     Ys = sort(Ys[:]);

     overlap = zeros(length(Ys)-1,1);
     Ns      = zeros(length(overlap),1);
     Sks     = zeros(length(overlap),1);
     Ms      = zeros(length(overlap),1);

     N = length(samples)/2;

     VarNull = 0;
     for i=1:length(overlap)
         sorted = sort([Ys[i],Ys[i+1]]);
         if sorted[1]>intervalMean[2]
         elseif sorted[2]<intervalMean[1]
         else
             overlap[i]=1;
             leftBounds = samplesSorted[:,2] .<= Ys[i];
             rightBounds = samplesSorted[:,1] .>= Ys[i+1];
             Ns[i] = sum(leftBounds) + sum(rightBounds);
             if Ns[i] == 0
                 VarNull = 1;
                 break;
             end
             leftOk   = samplesSorted[:,2] .* leftBounds;
             rightOk  = samplesSorted[:,1] .* rightBounds;

             #left = leftOk(leftOk ~=0);
             #right = rightOk(rightOk ~=0);
             left = filter!(x->x≠0,leftOk)
             right = filter!(x->x≠0,rightOk)
             Sks[i] = sum(left) + sum(right);
             test = Sks[i]/Ns[i];
             if sorted[1]<= test && sorted[2] >= test
                 if intervalMean[1]<= test && intervalMean[2] >= test
                     Ms[i] = (sum(left.^2) + sum(right.^2))/(length(samplesSorted)/2);
                 end
             end
         end
     end

     # Computing Lower Bound of Variance
     if VarNull == 0
         MsOks      = ones(length(Ms),1) .* Ms;
         SksFinal0  = Sks .* MsOks;
         NsFinal0   = Ns .* MsOks;

         #MsFinal1   = Ms(Ms ~= 0);
         #SksFinal1  = SksFinal0(SksFinal0 ~=0)./MsFinal1;
         #NsFinal1   = NsFinal0(NsFinal0 ~=0)./MsFinal1;

         MsFinal1 = filter!(x->x≠0,vec(Ms));
         SksFinal1  = filter!(x->x≠0,vec(SksFinal0))./MsFinal1;
         NsFinal1   = filter!(x->x≠0.0,vec(NsFinal0))./MsFinal1;

         vars = MsFinal1 - SksFinal1.^2 ./((N .* NsFinal1));
         intervalVariance[1] = minimum(vars) * (N/(N-1));
     end

     # Computing Upper Bound of Variance
     # Better way to cart prod julia?
     # p = num2cell(samples,2);
     # cart = cartprod(p{:});
     # a = collect(Iterators.product(skinny[1,:],skinny[2,:],skinny[3,:]))

     cart = cartesianProd(samples);

     vars = var.(cart');
     intervalVariance[2] = maximum(vars);

     return intervalMean , intervalVariance
end


# Cartesian Product algorithm
function cartesianProd(Sets1)

    results = [[]];

    SetNum = length(Sets1[:,1]);

    for i=1:SetNum
        currentSubArray  = Sets1[i,:];
        temp = [];
        for j=1:length(results)
            for k =1:length(currentSubArray)
                push!(temp, vcat(results[j],currentSubArray[k]));
            end
        end
        results = temp;
    end
    return results;
end


#=
kk = [interval(-1,0.5), interval(1,4),interval(2,4), interval(5,6)];

vl = LowerBoundVariance(kk)
=#
# True vals:
# Variance lower bound: 0.080619
# Variance upper bound: 3.55556
