using ProbabilityBoundsAnalysis, PyPlot

figSize = (8,8)
fontsize = 22

fig = figure("fig5_dist1", figsize = figSize)
x = U(0,3)
plot(x, name = "fig5_dist1", col = "red", fontsize = fontsize)
plt.axis("off")
savefig("fig5_dist1.pdf")

fig = figure("fig5_dist2",figsize = figSize)
y = N(0,1)
plot(y, name = "fig5_dist2", col = "blue", fontsize = fontsize)
plt.axis("off")
savefig("fig5_dist2.pdf")

fig = figure("fig5_pbox",figsize = figSize)
z = convFrechet(x,y, op =+)
plot(z, name = "fig5_pbox", col = "purple", fontsize = fontsize)
plt.axis("off")
savefig("fig5_pbox.pdf")
