###
#  This is a script for estimating descriptive statistics from interval
#  data sets, as outlined by S. Ferson (2007).
###



##
#   Has not yet been converted to the new format. May still work.
##


include("IntervalStatistics.jl");
using Plots;
using StatsBase;
using Distributions;

plotly();

# Skinny Data
skinny = [1.00 1.51; 2.68 2.98; 7.52 7.67; 7.52 8.35; 9.44 9.99; 3.66 4.58];
skinMean, skinVar = IntervalStatistics(skinny);

println(skinMean);
println(skinVar);

# Puffy Data
puffy = [3.5 6.4; 6.9 8.8; 6.1 8.4; 2.8 6.7; 3.5 9.7; 6.5 9.9; 0.15 3.8; 4.5 4.9; 7.1 7.9];
pufMean, pufVar = IntervalStatistics(puffy);


println(pufMean)
println(pufVar)

# Plotting skinny data
n = length(skinny)/2
plot([0,sort(skinny[:,1])], [(1:n)/n], linetype=:steppost);
plot!([0,sort(skinny[:,2])], [(1:n)/n], linetype=:steppost);

ran = skinMean[1]-2*sqrt(skinVar[2]):0.01:skinMean[2]+2*sqrt(skinVar[2]);

cdfupper = zeros(length(ran),1);
cdflower = zeros(length(ran),1);

upperDis1 = Normal(skinMean[1],sqrt(skinVar[1]));
upperDis2 = Normal(skinMean[1],sqrt(skinVar[2]));
lowerDis1 = Normal(skinMean[2],sqrt(skinVar[1]));
lowerDis2 = Normal(skinMean[2],sqrt(skinVar[2]));

plot!(ran, min.(cdf.(lowerDis1,ran),cdf.(lowerDis2,ran)));
plot!(ran, max.(cdf.(upperDis1,ran),cdf.(upperDis2,ran)));


 # Plotting puffy data
n = length(puffy)/2;
plot([0,sort(puffy[:,1])], [0,(1:n)/n], linetype=:steppost);
plot!([0,sort(puffy[:,2])], [0,(1:n)/n], linetype=:steppost);

ran = pufMean[1]-2*sqrt(pufVar[2]):0.01:pufMean[2]+2*sqrt(pufVar[2]);

cdfupper = zeros(length(ran),1);
cdflower = zeros(length(ran),1);

upperDis1 = Normal(pufMean[1],sqrt(pufVar[1]));
upperDis2 = Normal(pufMean[1],sqrt(pufVar[2]));
lowerDis1 = Normal(pufMean[2],sqrt(pufVar[1]));
lowerDis2 = Normal(pufMean[2],sqrt(pufVar[2]));

plot!(ran, min.(cdf.(lowerDis1,ran),cdf.(lowerDis2,ran)));
plot!(ran, max.(cdf.(upperDis1,ran),cdf.(upperDis2,ran)));
