%% This script generates the moments of the binary mixture using the linear noise approximation

addpath('../Common/');

syms kT v chi phi n1 s1 r V nA n alpha beta gamma k1 k2 c1 c2;

%Diffsion constant alpha subsumes factor 6

phi1 = n1*nA/(n1*nA + s1);
phi2 = (n - n1)*nA*v/(V - (n1*nA + s1)*v);
f = kT/v*(chi*phi*(1-phi)) + kT/v*(phi/nA*log(phi) + (1-phi)*log(1-phi));
f1 = simplify(subs(f, phi, phi1));
f2 = simplify(subs(f, phi, phi2));
%freeEnergy = simplify((n1*nA + s1)*v*f1 + (V - (n1*nA + s1)*v)*f2) + gamma*(V - (n1*nA + s1)*v)^(2/3);

df = simplify(diff(f, phi));
mu1 = simplify(subs(df, phi, phi1));
mu2 = simplify(subs(df, phi, phi2));


dmu = simplify(v*nA*(mu1 - mu2)); %Chemical potential difference
dLP = -2/3*v*36*pi^(1/3)*gamma/(V - (n1*nA + s1)*v)^(1/3); %Laplace pressure
dPi = simplify(v*(f1 - f2) - v*(phi1*mu1 - phi2*mu2) + dLP); %Osmotic pressure difference

w1 = simplify(alpha*n1/((n1*nA+s1)*v)^(2/3)*exp(-dmu/kT));
w2 = simplify(alpha*n1/((n1*nA+s1)*v)^(2/3));
%s2 = (V/v-nA*n*v)-s1;

w3 = simplify(beta*exp(-dPi/kT));%*(1/2+1/2*tanh(10*(s2-1)));
w4 = simplify(beta); %also 1/(s1+1) works

w5 = k1*r*(n1*nA + s1)*v; %syntehsis rate depends on volume area of the respective rate. Synthesis is also possible in the droplet.
w6 = k1*r*(V - (n1*nA + s1)*v);
w7 = k2*n1;
w8 = k2*(n-n1);

w9 = c1;
w10 = c2*r;

%     n1 s1   n  r
S = [ 1  -nA  0  0
    -1   nA   0  0
    0   1     0  0
    0  -1     0  0
    
    1  -nA    1  0
    0  0      1  0
    -1  nA    -1 0
    0   0     -1 0
    
    0   0     0  1
    0  0      0  -1];

params = {'kT' 'v' 'chi' 'V' 'nA' 'alpha' 'beta', 'gamma' 'k1', 'k2', 'c1', 'c2'};


h = { w1
    w2
    w3
    w4
    w5
    w6
    w7
    w8
    w9
    w10};


GenerateLNAODE(S, h, params, {'n1', 's1', 'n', 'r'}, 'ActiveBinaryMixture');

