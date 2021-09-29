using ProbabilityBoundsAnalysis, PyPlot, LaTeXStrings

figSize = (8,8)
fontsize = 24
fontsize_T = 36

x = uniform(1, interval(2,3))

fig = figure("fig8_pbox1", figsize = figSize)
plot(x, name = "fig8_pbox1", fontsize = fontsize)
PyPlot.title(L"$X \sim U(1, [2,3])$", fontsize = fontsize_T)
PyPlot.xlabel("X", fontsize = fontsize)
savefig("fig8_pbox1.pdf")

x2 = x^2

fig = figure("fig8_pbox2", figsize = figSize)
plot(x2, name = "fig8_pbox2", fontsize = fontsize)
PyPlot.title(L"$X^2$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig8_pbox2.pdf")

x3 = x^(0.5)

fig = figure("fig8_pbox3", figsize = figSize)
plot(x3, name = "fig8_pbox3", fontsize = fontsize)
PyPlot.title(L"$\sqrt{X}$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig8_pbox3.pdf")

x4 = exp(x)

fig = figure("fig8_pbox4", figsize = figSize)
plot(x4, name = "fig8_pbox4", fontsize = fontsize)
PyPlot.title(L"$e^{X}$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig8_pbox4.pdf")

x5 = log(x)

fig = figure("fig8_pbox5", figsize = figSize)
plot(x5, name = "fig8_pbox5", fontsize = fontsize)
PyPlot.title(L"$ln(X)$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig8_pbox5.pdf")

x6 = sin(x)

fig = figure("fig8_pbox6", figsize = figSize)
plot(x6, name = "fig8_pbox6", fontsize = fontsize)
PyPlot.title(L"$sin(X)$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig8_pbox6.pdf")

x7 = cos(x)

fig = figure("fig8_pbox7", figsize = figSize)
plot(x7, name = "fig8_pbox7", fontsize = fontsize)
PyPlot.title(L"$cos(X)$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig8_pbox7.pdf")


x8 = tanh(x)

fig = figure("fig8_pbox8", figsize = figSize)
plot(x8, name = "fig8_pbox8", fontsize = fontsize)
PyPlot.title(L"$tanh(X)$", fontsize = fontsize_T)
PyPlot.xlabel("Z", fontsize = fontsize)
savefig("fig8_pbox8.pdf")
