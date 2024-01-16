clear;
close all

addpath('../Common');
syms V nStar n1 n r kD k1 k2 c1 c2;

w1 = nStar*kD;
w2 = n1*kD;

w3 = k1*r*V;  %syntehsis rate depends on volume area of the respective rate.
w4 = k2*n1;
w5 = k2*(n-n1);

w6 = c1;
w7 = c2*r;

%     n1   n  r
S = [ 1    0  0
    -1     0  0
    
    1      1  0
    -1     -1 0
    0      -1 0
    
    0       0  1
    0       0  -1];

params = {'V', 'kD', 'nStar', 'k1', 'k2', 'c1', 'c2'};


h = { w1
    w2
    w3
    w4
    w5
    w6
    w7};


GenerateLNAODE(S, h, params, {'n1', 'n', 'r'}, 'ActiveBinaryMixture_SmallDropletLimit');