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

function plot(s ::pbox, fill = true; name = missing, col = missing, heading = missing, plotting = true, save = false, alpha = 0.2, fontsize = 12)

    if (isvacuous(s)); throw(ArgumentError("Pbox is vacuous"));end

    col1 = "red"; col2 = "black"; fillcol = "grey"
    if !(ismissing(col)); col1 = col2 = fillcol = col;end

    if !plotting; ioff();end
    if (ismissing(name)); fig = figure(figsize=(10,10)); else fig = figure(name,figsize=(10,10));end

    ax = fig.add_subplot()
    j = (0:(s.n-1))/s.n;

    PyPlot.step([s.u[:];s.u[s.n];s.d[s.n]], [j;1;1], color = col1, where = "pre");

    i = (1:(s.n))/s.n;
    PyPlot.step([s.u[1];s.d[1];s.d[:]], [0;0;i], color = col2,     where = "post");

    if fill
        Xs, Ylb, Yub = prepFillBounds(s);
        ax.fill_between(Xs, Ylb, Yub, alpha=alpha, color =fillcol)
    end
    if !(ismissing(heading)); title(heading, fontsize = fontsize);end
    xticks(fontsize = fontsize); yticks(fontsize = fontsize)
    xlabel("Distribution range",fontsize = fontsize); ylabel("CDF",fontsize=fontsize);

    if save; savefig("$name.png"); close(fig);end
    ion()
end

plot(s :: Union{<:Real, Interval{<:Real}}, fill = true; name = missing, col = missing, plotting = true, alpha = 0.2, fontsize = 12) = plot(makepbox(s), fill, name = name, col = col, plotting = plotting, alpha = alpha, fontsize = fontsize)


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

###
#       Copulas and bivariate p-box plots
###

function plotStep(cop :: copula; name = missing, pn = 50, fontsize = 18, alpha = 0.3)

    AU = cop.cdfU
    AD = cop.cdfD
    m = size(AU)[1];
    if m < pn; ppn = m; else ppn = pn; end      #pn is the plot "sub-sample", we can't plot all elements

    x = y = range(0, stop = 1,length=m)
    nm = round(m/ppn);

    x = x[1:Int(nm):end]
    oneIn = x[end] == 1

    if !oneIn; x = [x; 1]; end                             # We always want to include 1 in the plot

    zU = AU[1:Int(nm):end,1:Int(nm):end]
    zD = AD[1:Int(nm):end,1:Int(nm):end]

    if !oneIn;

        zU = [zU AU[end, 1:Int(nm):end]]
        zU = [zU' [AU[1:Int(nm):end, end]; 1]]

        zD = [zD AD[end, 1:Int(nm):end]]
        zD = [zD' [AD[1:Int(nm):end, end]; 1]]
    end

    if ismissing(name); fig = figure(figsize=(10,10)) else fig = figure(name,figsize=(10,10));end
    ax = fig.add_subplot(1,1,1,projection="3d")

    Ufaces = prepPolysU(zU,x);
    Dfaces = prepPolysD(zD,x);

    p3Ufaces = PyObject(art3d.Poly3DCollection(Ufaces, alpha=alpha))
    p3Dfaces = PyObject(art3d.Poly3DCollection(Dfaces, alpha=1))

    face_colorU = [1, 0, 0]
    face_colorD = [0, 0, 1]
    edge_color = [0.4, 0.4, 0.4]

    pycall(p3Ufaces.set_facecolor, PyAny, face_colorU)
    pycall(p3Ufaces.set_edgecolor, PyAny, edge_color)

    pycall(p3Dfaces.set_facecolor, PyAny, face_colorD)
    pycall(p3Dfaces.set_edgecolor, PyAny, edge_color)

    pycall(ax.add_collection3d, PyAny, p3Dfaces)
    pycall(ax.add_collection3d, PyAny, p3Ufaces)


end

function prepPolysU(us, x)

    ppn = length(x)

    p1 = [x[1]; x[1]; us[2, 2]]
    p2 = [x[1]; x[2]; us[2, 2]]
    p3 = [x[2]; x[2]; us[2, 2]]
    p4 = [x[2]; x[1]; us[2, 2]]
    tops = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]

    for i = 1:1:ppn-1
        for j = 1:1:ppn-1
            p1 = [x[i]; x[j]; us[i+1,j+1]]
            p2 = [x[i]; x[j+1]; us[i+1, j+1]]
            p3 = [x[i+1]; x[j+1]; us[i+1, j+1]]
            p4 = [x[i+1]; x[j]; us[i+1, j+1]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            tops = [tops; this]
        end
    end

    tops = tops[2:end]

    p1 = [x[1]; x[1]; us[1, 1]]
    p2 = [x[1]; x[1]; us[2, 2]]
    p3 = [x[1]; x[2]; us[2, 2]]
    p4 = [x[1]; x[2]; us[1, 1]]

    sides = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]

    for i = 1:1:ppn-1

        for j = i:1:ppn-1

            p1 = [x[i]; x[j]; us[i,j+1]]
            p2 = [x[i]; x[j]; us[i+1, j+1]]
            p3 = [x[i]; x[j+1]; us[i+1, j+1]]
            p4 = [x[i]; x[j+1]; us[i, j+1]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]

            p1 = [x[i]; x[j]; us[i+1,j]]
            p2 = [x[i]; x[j]; us[i+1, j+1]]
            p3 = [x[i+1]; x[j]; us[i+1, j+1]]
            p4 = [x[i+1]; x[j]; us[i+1, j]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]

            p1 = [x[j]; x[i]; us[j,i+1]]
            p2 = [x[j]; x[i]; us[j+1, i+1]]
            p3 = [x[j]; x[i+1]; us[j+1, i+1]]
            p4 = [x[j]; x[i+1]; us[j, i+1]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]

            p1 = [x[j]; x[i]; us[j+1,i]]
            p2 = [x[j]; x[i]; us[j+1, i+1]]
            p3 = [x[j+1]; x[i]; us[j+1, i+1]]
            p4 = [x[j+1]; x[i]; us[j+1, i]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]

        end
    end

    sides = unique(sides, dims=1)

    return [tops; sides]

end

function prepPolysD(us, x)

    ppn = length(x)

    p1 = [x[1]; x[1]; us[1, 1]]
    p2 = [x[1]; x[2]; us[1, 1]]
    p3 = [x[2]; x[2]; us[1, 1]]
    p4 = [x[2]; x[1]; us[1, 1]]
    tops = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]


    for i = 1:1:ppn-1
        for j = i:1:ppn-1
            p1 = [x[i]; x[j]; us[i+1, j]]
            p2 = [x[i]; x[j+1]; us[i+1, j]]
            p3 = [x[i+1]; x[j+1]; us[i+1, j]]
            p4 = [x[i+1]; x[j]; us[i+1, j]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            tops = [tops; this]

            p1 = [x[j]; x[i]; us[j, i+1]]
            p2 = [x[j]; x[i+1]; us[j, i+1]]
            p3 = [x[j+1]; x[i+1]; us[j, i+1]]
            p4 = [x[j+1]; x[i]; us[j, i+1]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            tops = [tops; this]

        end
    end

    p1 = [x[2]; x[2]; us[1, 1]]
    p2 = [x[2]; x[2]; us[2, 2]]
    p3 = [x[2]; x[3]; us[2, 2]]
    p4 = [x[2]; x[3]; us[1, 1]]

    sides = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]

    for i = 1:1:ppn-1
        for j = i:1:ppn-1

            p1 = [x[i]; x[j]; us[i+1,j]]
            p2 = [x[i]; x[j]; us[i, j]]
            p3 = [x[i]; x[j+1]; us[i, j]]
            p4 = [x[i]; x[j+1]; us[i+1, j]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]

            p1 = [x[j]; x[i]; us[j,i+1]]
            p2 = [x[j]; x[i]; us[j, i]]
            p3 = [x[j]; x[i+1]; us[j, i]]
            p4 = [x[j]; x[i+1]; us[j, i+1]]
            this = [[tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)]]
            sides = [sides; this]

        end
    end


    tops = unique(tops, dims=1)
    sides = unique(sides, dims=1)

    return [tops; sides]

end

function plot(x :: copula; name = missing, pn = 50, fontsize = 18, alpha = 0.7)

    AU = x.cdfU
    AD = x.cdfD
    m = size(AU)[1];
    if m < pn; ppn = m; else ppn = pn; end

    x = y = range(0, stop = 1,length=ppn)
    xgrid = repeat(x',ppn,1)
    ygrid = repeat(y,1,ppn)

    nm = round(m/ppn);

    zU = AU[1:Int(nm):end,1:Int(nm):end]
    zD = AD[1:Int(nm):end,1:Int(nm):end]

    if ismissing(name); fig = figure(figsize=(10,10)) else fig = figure(name,figsize=(10,10));end
    ax = fig.add_subplot(1,1,1,projection="3d")

    plot_surface(xgrid, ygrid, zD, rstride=2,edgecolors="b", cstride=2, alpha=1, linewidth=1, cmap=ColorMap("Blues"))
    plot_surface(xgrid, ygrid, zU, rstride=2,edgecolors="r", cstride=2, alpha=alpha, linewidth=1, cmap=ColorMap("RdGy"))

    ax.view_init(45-27, 180+ 26)

    xlabel("U",fontsize = fontsize)
    ylabel("V",fontsize = fontsize)
    ax.zaxis.set_rotate_label(false);
    zlabel("C(u,v)", rotation = 0, fontsize = fontsize)
    xticks(fontsize = fontsize ÷ 1.3); yticks(fontsize = fontsize ÷ 1.3);
    tight_layout()
end

function plot(x :: bivpbox; name = missing, pn = 50, fontsize = 18, alpha = 0.7)

    AU = x.C.cdfU
    AD = x.C.cdfD
    m = size(AU)[1];
    if m < pn; ppn = m; else ppn = pn; end

    xus = x.marg1.u; xds = x.marg1.d;
    yus = x.marg2.u; yds = x.marg2.d;

    nm = round(m/ppn);

    xus = xus[1:Int(nm):end]; xds = xds[1:Int(nm):end];
    yus = yus[1:Int(nm):end]; yds = yds[1:Int(nm):end];

    zU = AU[1:Int(nm):end,1:Int(nm):end]
    zD = AD[1:Int(nm):end,1:Int(nm):end]

    uGridx, uGridy = (repeat(xus',length(yus),1),repeat(yus,1,length(xus)))
    dGridx, dGridy = (repeat(xds',length(yds),1),repeat(yds,1,length(xds)))

    if ismissing(name); fig = figure(figsize=(10,10)) else fig = figure(name,figsize=(10,10));end
    ax = fig.add_subplot(1,1,1,projection="3d")


    plot_surface(dGridx, dGridy, zD, rstride=2,edgecolors="b", cstride=2, alpha=1, linewidth=1, cmap=ColorMap("Blues"))
    plot_surface(uGridx, uGridy, zU, rstride=2,edgecolors="r", cstride=2, alpha=alpha, linewidth=1, cmap=ColorMap("RdGy"))

    ax.view_init(45-27, 180+ 26)

    xlabel("y",fontsize = fontsize)
    ylabel("x",fontsize = fontsize)
    ax.zaxis.set_rotate_label(false);
    zlabel("H(x,y)", rotation = 0, fontsize = fontsize)
    tight_layout()
end

slice(x,y) = pycall(pybuiltin("slice"), PyObject, x, y)         # For performing python array slicing for

function scatter(a :: Array{Float64,2}; title = "samples")

    x = a[:,1];
    y = a[:,2];

    fig = plt.figure(title, figsize=(10, 10))
    grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
    main_ax = fig.add_subplot(get(grid, (slice(1,4),slice(0,3))))
    y_hist  = fig.add_subplot(get(grid, (slice(1,4),3)), xticklabels=[], sharey=main_ax)
    x_hist  = fig.add_subplot(get(grid, (0,slice(0,3))), yticklabels=[], sharex=main_ax)

    # scatter points on the main axes
    main_ax.plot(x, y, "o", markersize=3, alpha=0.2)

    # histogram on the attached axes
    x_hist.hist(x, 40, histtype="stepfilled", orientation="vertical", color="gray")

    y_hist.hist(y, 40, histtype="stepfilled", orientation="horizontal", color="gray")

end

function plotBoxes(xs :: Array{Interval{T},1}, ys  :: Array{Interval{T},1},  subpl = missing; linewidth = 1, alpha=0.2, fillcol= "grey") where T <: Real

    if ismissing(subpl);
        fig = plt.figure("boxes", figsize=(10, 10));
        subpl = fig.add_subplot()
    end
    Asize = length(xs)
    for i = 1:Asize

        xlo = xs[i].lo; xhi = xs[i].hi
        ylo = ys[i].lo; yhi = ys[i].hi

        subpl.plot([xlo; xhi], [ ylo; ylo], color = "red", linewidth = linewidth)
        subpl.plot([xhi; xhi], [ylo; yhi], color = "red", linewidth = linewidth)
        subpl.plot([xhi; xlo], [yhi; yhi], color = "red", linewidth = linewidth)
        subpl.plot([xlo; xlo], [yhi; ylo], color = "red", linewidth = linewidth)

        subpl.fill_between([xlo, xhi], [ylo, ylo], [yhi, yhi], alpha=alpha, color =fillcol)
    end

end

function scatter(a :: Array{Interval{T},2}; title = "samples", fontsize = 18) where T <: Real

    x = a[:,1];
    y = a[:,2];

    fig = plt.figure(title, figsize=(10, 10))
    grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
    main_ax = fig.add_subplot(get(grid, (slice(1,4),slice(0,3))))

    plotBoxes(x, y, main_ax,linewidth = 0.5, alpha=0.1, fillcol= "grey")

    xticks(fontsize = fontsize); yticks(fontsize = fontsize)
    xlabel("x",fontsize = fontsize); ylabel("y",fontsize=fontsize);

    y_hist  = fig.add_subplot(get(grid, (slice(1,4),3)), xticklabels=[], sharey=main_ax)
    x_hist  = fig.add_subplot(get(grid, (0,slice(0,3))), yticklabels=[], sharex=main_ax)

    xPbox = pbox(x);
    yPbox = pbox(y);

    alpha = 0.2;

    col1 = "red"; col2 = "black"; fillcol = "grey"

    j = (0:(xPbox.n-1))/xPbox.n;
    x_hist.step([xPbox.u[:];xPbox.u[xPbox.n];xPbox.d[xPbox.n]], [j;1;1], color = col1, where = "pre");

    i = (1:(xPbox.n))/xPbox.n;
    x_hist.step([xPbox.u[1];xPbox.d[1];xPbox.d[:]], [0;0;i], color = col2,     where = "post");

    Xs, Ylb, Yub = prepFillBounds(xPbox);
    x_hist.fill_between(Xs, Ylb, Yub, alpha=alpha, color =fillcol)

    j = (0:(yPbox.n-1))/yPbox.n;
    y_hist.step([j;1;1], [yPbox.u[:];yPbox.u[yPbox.n];yPbox.d[yPbox.n]], color = col1, where = "post");

    i = (1:(yPbox.n))/yPbox.n;
    y_hist.step([0;0;i],[yPbox.u[1];yPbox.d[1];yPbox.d[:]], color = col2,     where = "pre");

end
