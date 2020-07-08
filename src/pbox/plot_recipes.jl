######
# This file is part of the pba.jl package.
#
#   Definition of plotting recipes
#
#           University of Liverpool
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
#
#   Port of R code pba.r by Scott Ferson and Jason O'Rawe, Applied Biomathematics
#   Origional code available at: https://github.com/ScottFerson/pba.r
######

function plot(s ::pbox, fill = true; name = missing, col = missing, alpha = 0.2, fontsize = 12)

    if (isvacuous(s)); throw(ArgumentError("Pbox is vacuous"));end
    if (ismissing(name)) name = s.id; end

    col1 = "red"; col2 = "black"; fillcol = "grey"
    if !(ismissing(col)); col1 = col2 = fillcol = col;end

    fig = figure(name,figsize=(10,10))
    ax = fig.add_subplot()
    j = (0:(s.n-1))/s.n;

    PyPlot.step([s.u[:];s.u[s.n];s.d[s.n]], [j;1;1], color = col1, where = "pre");

    i = (1:(s.n))/s.n;
    PyPlot.step([s.u[1];s.d[1];s.d[:]], [0;0;i], color = col2,     where = "post");

    if fill
        Xs, Ylb, Yub = prepFillBounds(s);
        ax.fill_between(Xs, Ylb, Yub, alpha=alpha, color =fillcol)
    end

    xticks(fontsize = fontsize); yticks(fontsize = fontsize)
    xlabel("Distribution range",fontsize = fontsize); ylabel("CDF",fontsize=fontsize);

end

#plot(s ::pbox, fill = true; name = missing, col = missing, alpha = 0.2, fontsize = 12) = plot(makepbox(s), fill, name=name, col=col, alpha=alpha, fontsize=fontsize)

###
#   Prepares bounds for use with fill between
###
function prepFillBounds(x)

    d = x.d; u = x.u;

    is = range(0,stop =1 , length = x.n+1)
    di = is[2:end]; ui = is[1:end-1];

    Xs = sort([d; d; u; u]);
    nums = length(Xs)

    Ylb = zeros(nums,1); Yub = zeros(nums,1);

    for i = 1:2:nums
        indLb = findfirst(Xs[i] .<= d)
        indUb = findlast(Xs[i]  .>= u)

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
            if Xs[i] ∈ d
                Yub[i] = di[indUb]
            else
                Yub[i] = ui[indUb]
            end
        else
            Yub[i] = 0
            Yub[i+1]=0
        end

    end

    return Xs, Ylb[:], Yub[:]
end

#=
function plotpbox(s ::pbox; cumulative=ProbabilityBoundsAnalysis.cumulative, name = missing, col = missing)
    if (!isvacuous(s))
        if (ismissing(name)) name = s.name; end

        if cumulative ylab = "Cumulative probability" ; else ylab = "Exceedance probability";end

        j = (0:(s.n-1))/s.n;

        if (ismissing(col)) col = :red;end
        plot([s.u[:];s.u[s.n];s.d[s.n]], [j;1;1], colour = :red,linetype=:steppre, xlabel = "Distribution Range", ylabel = "CDF",legend = false,);

        i = (1:(s.n))/s.n;
        plot!([s.u[1];s.d[1];s.d[:]], [0;0;i], colour = :black,linetype=:steppost,legend = false);

        #plot!([s.u[:];s.u[s.n];s.d[s.n]], [j;1;1], fill = ([s.u[1];s.d[1];s.d[:]], [0;0;i], 0.5, :grey))


    end
end
=#



#=
@recipes function f(x :: AbstractPbox)


    seriestype := :steppre;
    arrow = :arrow;
    linealpha --> 0.5;
    linewidth = 4;


end

@recipes function f(x :: AbstractInterval)

end
=#
