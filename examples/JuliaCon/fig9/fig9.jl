using ProbabilityBoundsAnalysis, PyPlot, LaTeXStrings

fontsize = 24
fontsize_T = 36
figSize = (8,8)

#x1 = uniform(1, interval(2,3))
#x2 = uniform(1, interval(2,3))

x1 = uniform(interval(1,2), 4)
x2 = uniform(interval(1,2), 4)

fig = figure("fig9_pbox1", figsize = figSize)

plot(x1, name = "fig9_pbox1", fontsize = fontsize)
PyPlot.title(L"$X,Y \sim U([1,2], 4)$", fontsize = fontsize_T)
PyPlot.xlabel("X, Y", fontsize = fontsize)
savefig("fig9_pbox1.pdf")

z1 = convIndep(x1,x2, op = +)

fig = figure("fig9_pbox2", figsize = figSize)
plot(z1, name = "fig9_pbox2", fontsize = fontsize, col = "red")
PyPlot.title(L"$Z = X + Y$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig9_pbox2.pdf")


z2 = convPerfect(x1,x2, op = +)

fig = figure("fig9_pbox3", figsize = figSize)
plot(z2, name = "fig9_pbox3", fontsize = fontsize, col = "blue")
PyPlot.title(L"$Z = X + Y$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig9_pbox3.pdf")


C = GauCopula(-0.5)

z3 = sigma(x1,x2, op = *, C = C)

fig = figure("fig9_pbox4", figsize = figSize)
plot(z3, name = "fig9_pbox4", fontsize = fontsize, col = "green")
PyPlot.title(L"$Z = X \times Y$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig9_pbox4.pdf")

ax1 = gca()
XLIMS = ax1.get_xlim()

z4 = convOpposite(x1, x2 ,op = *)

fig = figure("fig9_pbox5", figsize = figSize)
plot(z4, name = "fig9_pbox5", fontsize = fontsize, col = "purple")
PyPlot.title(L"$Z = X \times Y$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
PyPlot.xlim(XLIMS)

savefig("fig9_pbox5.pdf")


z5 = convPerfect(x1, x2, op = /)

fig = figure("fig9_pbox6", figsize = figSize)
plot(z5, name = "fig9_pbox6", fontsize = fontsize, col = "blue")
PyPlot.title(L"$Z = X \div Y$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig9_pbox6.pdf")


C = GauCopula(0.5)

z6 = sigma(x1, x2, op = ^, C = C)

fig = figure("fig9_pbox7", figsize = figSize)
plot(z6, name = "fig9_pbox7", fontsize = fontsize, col = "orange")
PyPlot.title(L"$Z = X^{Y}$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig9_pbox8.pdf")


z7 = convIndep(x1, x2, op = /)

fig = figure("fig9_pbox8", figsize = figSize)
plot(z7, name = "fig9_pbox8", fontsize = fontsize, col = "red")
PyPlot.title(L"$Z = X \div Y}$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig9_pbox7.pdf")



#=
z2 = convIndep(x1,x2, op = *)

fig = figure("fig9_pbox3", figsize = figSize)
plot(z2, name = "fig9_pbox3", fontsize = fontsize, col = "red")
PyPlot.title(L"$Z = X \times Y$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig9_pbox3.pdf")
=#
