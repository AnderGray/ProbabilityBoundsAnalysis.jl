using ProbabilityBoundsAnalysis, PyPlot

a = normal(10,1);
b = uniform(1,2);

plot(a + b, fontsize = 22)
PyPlot.xlabel("a + b", fontsize = 26)
savefig("sumFrechet.png")
plot(a - b, fontsize = 22)
PyPlot.xlabel("a - b", fontsize = 26)
savefig("subFrechet.png")
plot(a * b, fontsize = 22)
PyPlot.xlabel("a * b", fontsize = 26)
savefig("prodFrechet.png")
plot(a / b, fontsize = 22)
PyPlot.xlabel("a / b", fontsize = 26)
savefig("divFrechet.png")

plot(convIndep(a, b), fontsize = 22)
PyPlot.xlabel("a + b (independent)", fontsize = 26)
savefig("sumIndep.png")
plot(convIndep(a, b, op = -), fontsize = 22)
PyPlot.xlabel("a - b (independent)", fontsize = 26)
savefig("subIndep.png")
plot(convIndep(a, b, op = *), fontsize = 22)
PyPlot.xlabel("a * b (independent)", fontsize = 26)
savefig("prodIndep.png")
plot(convIndep(a, b, op = /), fontsize = 22)
PyPlot.xlabel("a / b (independent)", fontsize = 26)
savefig("divIndep.png")

plot(convPerfect(a, b), fontsize = 22)
PyPlot.xlabel("a + b (perfect)", fontsize = 26)
savefig("sumPerfect.png")
plot(convPerfect(a, b, op = *), fontsize = 22)
PyPlot.xlabel("a * b (perfect)", fontsize = 26)
savefig("prodPerfect.png")
plot(convOpposite(a, b, op = +), fontsize = 22)
PyPlot.xlabel("a + b (opposite)", fontsize = 26)
savefig("sumOpposite.png")
plot(convOpposite(a, b, op = *), fontsize = 22)
PyPlot.xlabel("a * b (opposite)", fontsize = 26)
savefig("prodOpposite.png")


plot(conv(a, b, op = *, corr = -0.5), fontsize = 22)
PyPlot.xlabel("a * b (r = -0.5)", fontsize = 26)
savefig("prodCor1.png")

plot(conv(a, b, op = *, corr = 0.5), fontsize = 22)
PyPlot.xlabel("a * b (r = 0.5)", fontsize = 26)
savefig("prodCor2.png")

C1 = Clayton(5)
C2 = Frank(-2)

plot(convCorr(a, b, op = *, C = C1), fontsize = 22)
PyPlot.xlabel("a * b (C = Clayton(5))", fontsize = 26)
savefig("prodCop1.png")

plot(convCorr(a, b, op = *, C = C2), fontsize = 22)
PyPlot.xlabel("a * b (C = Frank(-2))", fontsize = 26)
savefig("prodCop2.png")
