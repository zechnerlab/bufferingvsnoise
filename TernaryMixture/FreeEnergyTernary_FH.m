function [F] = FreeEnergyTernary_FH(a1, b1, s1, chiAA, chiAB, chiBB, kT, v, nA, V, a, b, gamma)

%Volume fractions and droplet size
phiA1 = a1*nA/((a1*nA + b1*nA) + s1);
phiB1 = b1*nA/((a1*nA + b1*nA) + s1);
V2 = V-(a1*nA + b1*nA + s1)*v;
phiA2 = (a-a1)*nA*v / V2;
phiB2 = (b-b1)*nA*v / V2;

%Droplet surface area A
A = 36*pi^(1/3)*V2^(2/3);

F = ((a1*nA + b1*nA) + s1)*v*FreeEnergyDensityTernary_FH(phiA1, phiB1, chiAA, chiAB, chiBB, kT, v, nA) +...
    V2*FreeEnergyDensityTernary_FH(phiA2, phiB2, chiAA, chiAB, chiBB, kT, v, nA) + gamma*A;


end