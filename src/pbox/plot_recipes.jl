######
# This file is part of the pba.jl package.
#
# Definition of plotting recipes
#
#           University of Liverpool
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
#
#   Port of R code pba.r by Scott Ferson
#   Origional code available at:
#   Based on the paper:
######

function plotpbox(s ::pbox; cumulative=pba.cumulative, name = missing, col = missing)
    if (!isvacuous(s))
        if (ismissing(name)) name = s.name; end

        if cumulative ylab = "Cumulative probability" ; else ylab = "Exceedance probability";end

        j = (0:(s.n-1))/s.n;

        if (ismissing(col)) col = :red;end
        plot([s.u[:];s.u[s.n];s.d[s.n]], [j;1;1], colour = :red,linetype=:steppre,xlabel = "Distribution Range", ylabel = "CDF",legend = false,);

        i = (1:(s.n))/s.n;
        plot!([s.u[1];s.d[1];s.d[:]], [0;0;i], colour = :black,linetype=:steppost,legend = false);

        #plot!([s.u[:];s.u[s.n];s.d[s.n]], [j;1;1], fill = ([s.u[1];s.d[1];s.d[:]], [0;0;i], 0.5, :grey))


    end
end



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
