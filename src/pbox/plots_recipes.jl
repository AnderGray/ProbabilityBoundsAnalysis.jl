######
# This file is part of the pba.jl package.
#
#   Definition of plotting recipes
#
#           University of Liverpool
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
#
#   About 50% of this file is a port of R code pba.r by Scott Ferson and Jason O'Rawe, Applied Biomathematics
#   Origional code available at: https://github.com/ScottFerson/pba.r
######


DEFAULT_ALPHA = 0.2
DEFAULT_LABEL = ""
DEFAULT_GRID = false
DEFAULT_LEGEND = false
DEFAULT_FILL = :gray
DEFAULT_COLOUR_UPPER = :red
DEFAULT_COLOUR_LOWER = :black

DEFAULT_FONT_SIZE = 18
DEFAULT_TICK_SIZE = 12


#fill = true; name = missing, col = missing, heading = missing, plotting = true, save = false, alpha = 0.2, fontsize = 12

@recipe function _pbox_plot(s ::pbox)

    j = (0:(s.n-1))/s.n;
    i = (1:(s.n))/s.n;
    
    grid --> DEFAULT_GRID
    legend --> DEFAULT_LEGEND
    ylabel --> "cdf"

    xguidefontsize --> DEFAULT_FONT_SIZE
    yguidefontsize --> DEFAULT_FONT_SIZE
    
    Xs, Ylb, Yub = prepFillBounds(s)

    @series begin
        fillrange := Yub
        color := DEFAULT_FILL
        fillalpha := DEFAULT_ALPHA
        Xs, Ylb
    end

    @series begin
        linetype := :steppre
        color --> DEFAULT_COLOUR_UPPER
        alpha := 1
        label := ""
        [s.u[:];s.u[s.n];s.d[s.n]], [j;1;1]
    end

    @series begin
        linetype := :steppost
        color --> DEFAULT_COLOUR_LOWER
        alpha := 1
        label := ""
        [s.u[1];s.d[1];s.d[:]], [0;0;i]
    end    
end

@recipe function _plot_interval(s :: Interval{<:Real})
    xlims --> (s.lo - diam(s)/5, s.hi + diam(s)/5)
    makepbox(s)
end

# recipes --> for converting data into other atributes, e.g. converting a bar series into a shape series
#

function prepFillBounds(x)

    d = x.d; u = x.u;

    is = range(0,stop =1 , length = x.n+1)
    di = is[2:end]; ui = is[1:end-1];

    Xs = sort([d; d; u; u]);
    nums = length(Xs)

    Ylb = zeros(nums,1); Yub = zeros(nums,1);

    for i = 1:2:nums

        indUb = findlast(Xs[i]  .>= u)
        indLb = findfirst(Xs[i] .<= d)

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
            if Xs[i] ∈ u
                Yub[i] = ui[indUb]
            else
                Yub[i] = di[indUb]
            end
        else
            Yub[i] = 0
            Yub[i+1]=0
        end

    end

    return Xs, Ylb[:], Yub[:]
end
