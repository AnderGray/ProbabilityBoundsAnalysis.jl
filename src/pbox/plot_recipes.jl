######
# This file is part of the pba.jl package.
#
# Definition of plotting recipes
#
#           University of Liverpool
######

function plotpbox(s ::pbox; cumulative=pba.cumulative, name = missing, col = missing)
    if (!isvacuous(s))
        if (ismissing(name)) name = s.name; end

        if cumulative ylab = "Cumulative probability" ; else ylab = "Exceedance probability";end

        i = (0:(s.n-1))/s.n;

        if (ismissing(col)) col = :red;end
        plot([s.u[:];s.u[s.n];s.d[s.n]], [i;1;1], colour = :red,linetype=:steppre);

        i = (1:(s.n))/s.n;
        plot!([s.u[1];s.d[1];s.d[:]], [0;0;i], colour = :black,linetype=:steppost);

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
