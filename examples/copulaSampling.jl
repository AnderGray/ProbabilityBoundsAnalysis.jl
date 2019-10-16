######
# This file is part of the pba.jl package.
#
#   Examples of how to sample from copulas
#
#           University of Liverpool
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
######

include("../src/pbox/copula.jl")

n = 200;


cop = π()           # Independent copula
plotCdf(cop);
wait_for_key();
sample(cop,10,plot=true)
PyPlot.clf()

cop = M()
plotCdf(cop);
wait_for_key()
sample(cop,2,plot=true)
PyPlot.clf()

cop = W()
plotCdf(cop);
wait_for_key()
sample(cop,2,plot=true)
PyPlot.clf()

cop = Gaussian(0)
plotCdf(cop);
wait_for_key()
PyPlot.clf()

cop = Gaussian(-0.8)
plotCdf(cop);
wait_for_key()
sample(cop,10,plot=true)
PyPlot.clf()

cop = Gaussian(0.8)
plotCdf(cop);
wait_for_key()
PyPlot.clf()

samps = sample(cop,1000);
scatter(samps[:,1], samps[:,2])
wait_for_key()
PyPlot.clf()

marginal1 = LogNormal(5,1);
marginal2 = Beta(10,2);

s1 = quantile.(marginal1, samps[:,1])
s2 = quantile.(marginal2, samps[:,2])

scatter(samps[:,1], samps[:,2])
wait_for_key()
PyPlot.clf()
