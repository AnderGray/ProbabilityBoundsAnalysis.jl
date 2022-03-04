######
# This file is part of the pba.jl package.
#
#   Examples of how to sample from copulas
#
#           University of Liverpool
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
######

using ProbabilityBoundsAnalysis, PyPlot

n = 200;


cop = πCop()           # Independent copula
plot(cop);
samps = sample(cop, n)
scatter(samps)
PyPlot.clf()

cop = M()
plot(cop);
samps = sample(cop, n)
scatter(samps)
PyPlot.clf()

cop = W()
plot(cop);
samps = sample(cop, n)
scatter(samps)
PyPlot.clf()

cop = Gaussian(0)
plot(cop);
samps = sample(cop, n)
scatter(samps)
PyPlot.clf()

cop = Gaussian(-0.8)
plot(cop);
samps = sample(cop, n)
scatter(samps)
PyPlot.clf()

cop = Gaussian(0.8)
plot(cop);
samps = sample(cop, n)
scatter(samps)
PyPlot.clf()

marginal1 = lognormal(5,1);
marginal2 = beta(10,2);

s1 = cut.(marginal1, samps[:,1])
s2 = cut.(marginal2, samps[:,2])

scatter(s1, s2)
PyPlot.clf()
