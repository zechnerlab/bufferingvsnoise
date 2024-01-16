function [F] = FreeEnergyBinary_FH_fun(n1, s1, chi, kT, v, nA, V, n, gamma)

%Volume fractions and droplet size
phi1 = n1*nA./(n1*nA + s1);
V2 = V-(n1*nA + s1)*v;
phi2 = (n-n1)*nA*v ./ V2;

%Droplet surface area A
A = 36*pi^(1/3)*V2.^(2/3);

F = (n1*nA + s1)*v.*FreeEnergyDensityBinary_FH_fun(phi1, chi, kT, v, nA) +...
    (V-(n1*nA + s1)*v).*FreeEnergyDensityBinary_FH_fun(phi2, chi, kT, v, nA) + gamma*A;

end