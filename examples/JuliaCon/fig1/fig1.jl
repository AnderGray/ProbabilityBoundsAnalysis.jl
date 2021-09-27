using ProbabilityBoundsAnalysis, PyPlot

ProbabilityBoundsAnalysis.setSteps(500)

a = normal(0, 1)                # mean = 0, var = 1
b = interval(-1,1)               # from IntervalArithmetic.jl
c = normal(interval(-1,1), 1)    # mean = [-1, 1] var = 1

fontsize = 24

fig = figure("fig1",figsize=(8,8))
plot(a, fontsize = fontsize, name = "fig1")
PyPlot.xlabel("X", fontsize = fontsize)
savefig("fig1_dist.pdf")

fig = figure("fig2",figsize=(8,8))
plot(b, fontsize = fontsize, name = "fig2")
PyPlot.xlabel("X", fontsize = fontsize)
PyPlot.xlim([-2, 2])
savefig("fig1_interval.pdf")

fig = figure("fig3",figsize=(8,8))
plot(c, fontsize = fontsize, name = "fig3")
PyPlot.xlabel("X", fontsize = fontsize)
savefig("fig1_pbox.pdf")
