####
#
# Frechet bounds with sampling
#                   University of Liverpool
###


using Distributions
using Plots
using StatsBase

function sample(Fx ::Sampleable{Univariate}, Fy ::Sampleable{Univariate}; N :: Int64= 100000, op = +, bining = collect(-10:0.01:10);)

#N = 100000  # number of samples

    u = rand(N)
    v = rand(N)

    #F = Uniform(0,1);
    #G = Uniform(0,1);
    #F = Normal(0,1);
    #G = Normal(0,1);

    x = quantile.(Fx,u);
    y = quantile.(Fy,v);

    if (op) ∈ (+,-)
        z = map(op,x,y);
    elseif (op) ∈ (*,/)
        z = map(op,x,y);
    end

    z = round.(z, digits=2);

    if (op) ∈ (+,*)
        Tau = max.(u .+ v .- 1, 0);
        Row = min.(u .+ v, 1);
    elseif op == -
        Tau = max.(u .- cdf.(Fy,-y), 0);
        Row = min.(u .- cdf.(Fy,-y), 1);
    elseif op == /
        Tau = max.(u .- cdf.(Fy,1 ./y),0);
        Row = min.(u .- cdf.(Fy,1 ./y), 1);
    end
    #Row =(u .+ v .- Tau);

    T = zeros(length(bining)-1);
    p = zeros(length(bining)-1);

    #row = inf(F(x) + G(y) - max(F(x)+G(y)-1,0))

    if (op) ∈ (+,*)

        for i = 1:length(bining)-1
            if !(isempty(Tau[z .== bining[i+1]]))
                T[i] = maximum(Tau[z .== bining[i+1]])
            end
            if !(isempty(Row[z .== bining[i]]))
                p[i] = minimum(Row[z .== bining[i]]);
            end
        end

    elseif (op) ∈ (-,/)

        for i = 1:length(bining)-1
            if !(isempty(Tau[z .== bining[i+1]]))
                T[i] = maximum(Tau[z .== bining[i+1]])
            end
            if !(isempty(Row[z .== bining[i]]))
                p[i] = 1 + minimum(Row[z .== bining[i]]);
            end
        end
    end

    plot(bining[2:end], T, colour = :black,linetype=:steppost);
    plot!(bining[1:end-1], p, colour = :red,linetype=:steppre);

end
