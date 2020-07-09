####
#   Examples of the tau rho convolutions
####

#include("../../src/pbaDebug.jl")

x = U(0,1);
y = U(0,1);

Cors = range(-1,stop = 1, length = 6)
colors = ["red", "blue", "purple", "black", "green", "orange"]

op = +;

for i = 1:length(Cors)

    Cop = GauCopula(Cors[i])

    if Cors[i] == -1; Cop = W();
elseif Cors[i] == 1; Cop = M(); 
else Cop = GauCopula(Cors[i]); end

    
    Fret1 = tauRho(x, y,Cop, op);

    plot(Fret1, name = "same", col = colors[i], alpha = 0.1)
    if i == length(Cors); savefig("SumUnif.png");end
end

