include("../src/pba.jl");
using Main.pba;


a = normal(interval(0,1),1);
b = normal(interval(0,1),1);

c = a + b;
