######
# This file is part of the pba.jl package.
#
#   Numerical experiments (Monte Carlo) for dependence between random variables
#   which have previously been used in an opertaion
#
#           University of Liverpool
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
#
#   For hints as to why this happens view section covariance alegra in: https://en.wikipedia.org/wiki/Algebra_of_random_variables
######

include("../src/pbox/dependency.jl")

#   Define Marginals and copula
a = π(Normal(10,2),Normal());

samps = sample(a,10^6);

x = samps[:,1]
y = samps[:,2]


# Sum
z = x .+y;

scatter(hcat(x,z), title = "Scatter x z (x+y)")
scatter(hcat(y,z), title = "Scatter y z (x+y)")

corZY = cor(z,y);
corZX = cor(z,x);

println("-----------------------------")
println("(x+y) correlation Z X: $corZX")
println("(x+y) correlation Z Y: $corZY")
println()

# Difference
z = x .-y;

scatter(hcat(x,z), title = "Scatter x z (x-y)")
scatter(hcat(y,z), title = "Scatter y z (x-y)")

corZY = cor(z,y);
corZX = cor(z,x);

println("-----------------------------")
println("(x-y) correlation Z X: $corZX")
println("(x-y) correlation Z Y: $corZY")
println()

# Product
z = x .*y;

scatter(hcat(x,z), title = "Scatter x z (x*y)")
scatter(hcat(y,z), title = "Scatter y z (x*y)")

corZY = cor(z,y);
corZX = cor(z,x);

println("-----------------------------")
println("(x*y) correlation Z X: $corZX")
println("(x*y) correlation Z Y: $corZY")
println()

# Quotient
z = x ./y;

scatter(hcat(x,z), title = "Scatter x z (x/y)")
scatter(hcat(y,z), title = "Scatter y z (x/y)")

corZY = cor(z,y);
corZX = cor(z,x);

println("-----------------------------")
println("(x/y) correlation Z X: $corZX")
println("(x/y) correlation Z Y: $corZY")
println()

# Inputs
corXY = cor(x,y)
scatter(hcat(x,y), title = "Scatter x y")
println("-----------------------------")
println("correlation X Y: $corXY")
println()
#=

# For repeated sums

N = 1;
z = x .+ y;

for i = 1:N

    global z = z .+ x;

    corZX = cor(z,x);

    #println("-----------------------------")
    #println("correlation after $i iterations: $corZX")
    #println()

end

scatter(hcat(x,z), title = "After $N iterations")

#scatter(hcat(x,z), title = "Scatter x z (x+y)")
#scatter(hcat(y,z), title = "Scatter y z (x+y)")

#corZY = cor(z,y);
#corZX = cor(z,x);
=#
