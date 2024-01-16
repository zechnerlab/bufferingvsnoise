%% Generates moments for the ternary mixture using the linear noise approximation

clear;
close all;

addpath('../Common/');

syms kT v chiAA chiAB chiBB phiA phiB a1 b1 s1 r V nA a b alpha beta gamma k1A k1B k2A k2B c1 c2;

%Diffsion constant alpha subsumes factor 6

phiA1 = a1*nA/(a1*nA + b1*nA + s1);
phiA2 = (a - a1)*nA*v/(V - (a1*nA + b1*nA + s1)*v);
phiB1 = b1*nA/(a1*nA + b1*nA + s1);
phiB2 = (b - b1)*nA*v/(V - (a1*nA + b1*nA + s1)*v);
f = kT/v*(chiAA*phiA*(1-phiA-phiB) + chiAB*phiA*phiB + chiBB*phiB*(1-phiA-phiB))...
    + kT/v*(phiA/nA*log(phiA) + phiB/nA*log(phiB) + (1-phiA-phiB)*log(1-phiA-phiB));
f1 = simplify(subs(f, [phiA, phiB], [phiA1, phiB1]));
f2 = simplify(subs(f, [phiA, phiB], [phiA2, phiB2]));

dfA = simplify(diff(f, phiA));
dfB = simplify(diff(f, phiB));
muA1 = simplify(subs(dfA, [phiA, phiB], [phiA1, phiB1]));
muA2 = simplify(subs(dfA, [phiA, phiB], [phiA2, phiB2]));
muB1 = simplify(subs(dfB, [phiA, phiB], [phiA1, phiB1]));
muB2 = simplify(subs(dfB, [phiA, phiB], [phiA2, phiB2]));

dmuA = simplify(v*nA*(muA1 - muA2)); %Chemical potential difference
dmuB = simplify(v*nA*(muB1 - muB2)); %Chemical potential difference
dLP = -2/3*v*36*pi^(1/3)*gamma/(V - (a1*nA + b1*nA + s1)*v)^(1/3); %Laplace pressure
dPi = simplify(v*(f1 - f2) - v*(phiA1*muA1 + phiB1*muB1 - phiA2*muA2 - phiB2*muB2) + dLP); %Osmotic pressure difference

w1 = simplify(alpha*a1/((a1*nA+b1*nA+s1)*v)^(2/3)*exp(-dmuA/kT));
w2 = simplify(alpha*a1/((a1*nA+b1*nA+s1)*v)^(2/3));

w3 = simplify(alpha*b1/((a1*nA+b1*nA+s1)*v)^(2/3)*exp(-dmuB/kT));
w4 = simplify(alpha*b1/((a1*nA+b1*nA+s1)*v)^(2/3));

w5 = simplify(beta*exp(-dPi/kT));%*(1/2+1/2*tanh(10*(s2-1)));
w6 = simplify(beta); %also 1/(s1+1) works

w7 = k1A*r*(a1*nA + b1*nA + s1)*v; %syntehsis rate depends on volume area of the respective rate. Synthesis is also possible in the droplet.
w8 = k1A*r*(V - (a1*nA + b1*nA + s1)*v);
w9 = k2A*a1;
w10 = k2A*(a-a1);

w11 = k1B*r*(a1*nA + b1*nA + s1)*v; %syntehsis rate depends on volume area of the respective rate. Synthesis is also possible in the droplet.
w12 = k1B*r*(V - (a1*nA + b1*nA + s1)*v);
w13 = k2B*b1;
w14 = k2B*(b-b1);

w15 = c1;
w16 = c2*r;


%     a1 b1   s1  a b  r
S = [ 1  0   -nA  0 0  0  
    -1  0   nA    0 0  0  
    0  1  -nA     0 0  0  
    0  -1  nA     0 0  0  
    
    0   0   1     0 0  0  
    0   0   -1    0 0  0  
    
    1  0 -nA    1 0 0  
    0  0 0      1 0 0  
    -1 0 nA    -1 0 0  
    0  0 0     -1 0 0  
    
    0  1 -nA    0  1 0 
    0  0 0      0  1 0 
    0 -1 nA     0 -1 0 
    0  0 0      0 -1 0 
    
    0  0 0      0 0  1  
    0  0 0      0 0 -1  
];

params = {'kT' 'v' 'chiAA' 'chiAB' 'chiBB' 'V' 'nA' 'alpha' 'beta', 'gamma' 'k1A', 'k1B', 'k2A', 'k2B', 'c1', 'c2'};


h = { w1
    w2
    w3
    w4
    w5
    w6
    w7
    w8
    w9
    w10
    w11
    w12
    w13
    w14
    w15
    w16
};


GenerateLNAODE(S, h, params, {'a1', 'b1', 's1', 'a', 'b', 'r'}, 'ActiveTernaryMixture_Correlation');
