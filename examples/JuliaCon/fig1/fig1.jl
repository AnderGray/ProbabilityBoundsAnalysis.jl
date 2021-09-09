using ProbabilityBoundsAnalysis, PyPlot

ProbabilityBoundsAnalysis.setSteps(500)

a = normal(0, 1)                # mean = 0, var = 1
b = interval(-1,1)               # from IntervalArithmetic.jl
c = normal(interval(-1,1), 1)    # mean = [-1, 1] var = 1

fontsize = 24

plot(a, fontsize = fontsize)
PyPlot.xlabel("X", fontsize = fontsize)
savefig("fig1_dist.pdf")
plot(b, fontsize = fontsize)
PyPlot.xlabel("X", fontsize = fontsize)
PyPlot.xlim([-2, 2])
savefig("fig1_interval.pdf")
plot(c, fontsize = fontsize)
PyPlot.xlabel("X", fontsize = fontsize)
savefig("fig1_pbox.pdf")
