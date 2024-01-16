function dN = ActiveBinaryMixture_SmallDropletLimit(t, N, V, kD, nStar, k1, k2, c1, c2)
n1=N(1);
n=N(2);
r=N(3);
C = reshape(N(4:end), 3, 3);

dX(1) = kD*nStar - kD*n1 - k2*n1 + V*k1*r;
dX(2) = V*k1*r - k2*n;
dX(3) = c1 - c2*r;

A(1, 1) = - k2 - kD;
A(1, 2) = 0;
A(1, 3) = V*k1;
A(2, 1) = 0;
A(2, 2) = -k2;
A(2, 3) = V*k1;
A(3, 1) = 0;
A(3, 2) = 0;
A(3, 3) = -c2;

BBT(1, 1) = k2*n1 + kD*n1 + kD*nStar + V*k1*r;
BBT(1, 2) = k2*n1 + V*k1*r;
BBT(1, 3) = 0;
BBT(2, 1) = k2*n1 + V*k1*r;
BBT(2, 2) = k2*n + V*k1*r;
BBT(2, 3) = 0;
BBT(3, 1) = 0;
BBT(3, 2) = 0;
BBT(3, 3) = c1 + c2*r;

dC = A*C + C*A' + BBT;
dN = [dX(:); reshape(dC, 9, 1)];
