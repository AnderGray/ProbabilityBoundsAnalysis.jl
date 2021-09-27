using ProbabilityBoundsAnalysis, PyPlot

ProbabilityBoundsAnalysis.setSteps(1000)

using3D()

x1 = beta(interval(4,8), 3)
x2 = normal(interval(0,2),2)

C = GauCopula(-0.8)

j1 = C(x1,x2)

ft = 24
ft_ticks = 20

name = "biv_cdf"
fig = figure(name,figsize=(10,10))
plot(j1, name = name, fontsize = ft)
PyPlot.zlabel("F(x,y)", fontsize = ft)

ax = gca()
PyPlot.xticks(fontsize = ft_ticks)
PyPlot.yticks(fontsize = ft_ticks)
#ax.zticks(fontsize = ft_ticks)
for t in ax.zaxis.get_major_ticks()
    t.label.set_fontsize(ft_ticks)
end

ax.view_init(30, -60)
fig.canvas.draw()
savefig("$name.pdf")

s = sample(j1, 100)
plotBoxes(s[:,1], s[:,2], alpha= 0.1)
PyPlot.xticks(fontsize = ft_ticks)
PyPlot.yticks(fontsize = ft_ticks)
PyPlot.xlabel("X", fontsize = ft)
PyPlot.ylabel("Y", fontsize = ft)

savefig("samples.pdf")
