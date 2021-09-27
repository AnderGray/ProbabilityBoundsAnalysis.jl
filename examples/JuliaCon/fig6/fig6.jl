using ProbabilityBoundsAnalysis, PyPlot, LaTeXStrings

figSize = (8,8)
fontsize = 22

x1 = MeanVar(0,1)
fig = figure("fig6_pbox1", figsize = figSize)
plot(x1, name = "fig6_pbox1", fontsize = fontsize)
PyPlot.title("MeanVar(0, 1)", fontsize = fontsize)
PyPlot.xlabel("X", fontsize = fontsize)
savefig("fig6_pbox1.pdf")

x2 = MeanMin(1,0)
fig = figure("fig6_pbox2", figsize = figSize)
plot(x2, name = "fig6_pbox2", fontsize = fontsize)
PyPlot.title("MeanMin(1, 0)", fontsize = fontsize)
PyPlot.xlabel("X", fontsize = fontsize)
savefig("fig6_pbox2.pdf")

x3 = MeanMinMax(1,0,2)
fig = figure("fig6_pbox3", figsize = figSize)
plot(x3, name = "fig6_pbox3", fontsize = fontsize)
PyPlot.title("MeanMinMax(1, 0, 2)", fontsize = fontsize)
PyPlot.xlabel("X", fontsize = fontsize)
savefig("fig6_pbox3.pdf")

x4 = MinMaxMeanVar(0,2,1,0.5)
fig = figure("fig6_pbox4", figsize = figSize)
plot(x4, name = "fig6_pbox4", fontsize = fontsize)
PyPlot.title("MinMaxMeanVar(0, 2, 1, 0.5)", fontsize = fontsize)
PyPlot.xlabel("X", fontsize = fontsize)
savefig("fig6_pbox4.pdf")
