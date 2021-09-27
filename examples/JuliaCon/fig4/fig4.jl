using ProbabilityBoundsAnalysis, PyPlot

figSize = (8,8)
fontsize = 22

fig = figure("fig4_dist", figsize = figSize)
x = exponential(0.6)
plot(x, name = "fig4_dist", col = "red", fontsize = fontsize)
plt.axis("off")
savefig("fig4_dist.pdf")

fig = figure("fig4_in",figsize = figSize)
y = interval(10, 20)
plot(y, name = "fig4_in", col = "blue", fontsize = fontsize)
PyPlot.xlim((8,22))
plt.axis("off")
savefig("fig4_in.pdf")

fig = figure("fig4_pbox",figsize = figSize)
z = x * makepbox(y)
plot(z, name = "fig4_pbox", col = "purple", fontsize = fontsize)
plt.axis("off")
savefig("fig4_pbox.pdf")
