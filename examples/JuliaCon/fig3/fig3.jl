using ProbabilityBoundsAnalysis, PyPlot, LaTeXStrings

ft = 24
ft2 = 30

a = normal(interval(0, 1), interval(2, 3))

b = uniform(interval(0, 1), 3)

c = beta(interval(0.7, 1), interval(3, 4))

d = gamma(interval(5, 6), 2)

e = lognormal(interval(2, 3), interval(1, 5))

f = exponential(interval(0.4, 0.6))

g = chisq(interval(20, 50))

h = cauchy(interval(1, 100), 1)


plot(a, fontsize = ft)
PyPlot.title(L"$X \sim N([0,1], [2,3])$", fontsize = ft2)
PyPlot.xlabel("X", fontsize = ft)
savefig("normal.pdf")

plot(b, fontsize = ft)
PyPlot.title(L"$X \sim U([0,1], 3)$", fontsize = ft2)
PyPlot.xlabel("X", fontsize = ft)
savefig("uniform.pdf")

plot(c, fontsize = ft)
PyPlot.title(L"$X \sim Beta([0.7,1], [3,4])$", fontsize = ft2)
PyPlot.xlabel("X", fontsize = ft)
savefig("beta.pdf")

plot(d, fontsize = ft)
PyPlot.title(L"$X \sim Gamma([5,6], 2)$", fontsize = ft2)
PyPlot.xlabel("X", fontsize = ft)
savefig("gamma.pdf")

plot(e, fontsize = ft)
PyPlot.title(L"$X \sim Lognormal([2,3], [1,5])$", fontsize = ft2)
PyPlot.xlabel("X", fontsize = ft)
savefig("lognormal.pdf")

plot(f, fontsize = ft)
PyPlot.title(L"$X \sim Exp([0.4,0.6])$", fontsize = ft2)
PyPlot.xlabel("X", fontsize = ft)
savefig("exponential.pdf")

plot(g, fontsize = ft)
PyPlot.title(L"$X \sim Chisquared([20,50])$", fontsize = ft2)
PyPlot.xlabel("X", fontsize = ft)
savefig("chisquared.pdf")

plot(h, fontsize = ft)
PyPlot.title(L"$X \sim Chauchy([1,100], 1)$", fontsize = ft2)
PyPlot.xlabel("X", fontsize = ft)
savefig("cauchy.pdf")
