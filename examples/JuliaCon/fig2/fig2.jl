using ProbabilityBoundsAnalysis, PyPlot, Distributions

fontsize = 24

xs = range(-3, 3, length = 500)
cdfsX = cdf.(Normal(0,1), xs)

ProbabilityBoundsAnalysis.setSteps(25)
a = normal(0, 1)                # mean = 0, var = 1

fig = figure("fig1",figsize=(8,8))
plot(xs,cdfsX, color = "blue", linewidth = 2)
plot(a, fontsize = fontsize, name = "fig1")
PyPlot.xlabel("X", fontsize = fontsize)
savefig("fig2_pbox1.pdf")

ProbabilityBoundsAnalysis.setSteps(5)
a = normal(0, 1)                # mean = 0, var = 1

fig = figure("fig2",figsize=(8,8))
plot(xs,cdfsX, color = "blue", linewidth = 2)
plot(a, fontsize = fontsize, name = "fig2")
PyPlot.xlabel("X", fontsize = fontsize)
savefig("fig2_pbox2.pdf")

ProbabilityBoundsAnalysis.setSteps(10)
a = normal(interval(-1,1), 1)                # mean = 0, var = 1

fig = figure("fig3",figsize=(8,8))
plot(xs .- 1,cdfsX,  color = "blue", linewidth = 2)
plot(xs .+ 1,cdfsX, color = "blue", linewidth = 2)
plot(a, fontsize = fontsize, name = "fig3")
PyPlot.xlabel("X", fontsize = fontsize)
savefig("fig2_pbox3.pdf")
