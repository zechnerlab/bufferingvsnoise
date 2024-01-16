function dN = ActiveBinaryMixture(t, N, kT, v, chi, V, nA, alpha, beta, gamma, k1, k2, c1, c2)
n1=N(1);
s1=N(2);
n=N(3);
r=N(4);
C = reshape(N(5:end), 4, 4);

dX(1) = (alpha*n1*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT))/(v*(s1 + n1*nA))^(2/3) - (alpha*n1)/(v*(s1 + n1*nA))^(2/3) - k2*n1 + k1*r*v*(s1 + n1*nA);
dX(2) = beta*exp((v*((kT*(n - n1)*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(V - v*(s1 + n1*nA)) - (kT*n1*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(v*(s1 + n1*nA))) - v*((kT*(log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA)))*((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1) - (v*log((nA*v*(n - n1))/(V - v*(s1 + n1*nA)))*(n - n1))/(V - v*(s1 + n1*nA))))/v + (kT*(s1*log(s1/(s1 + n1*nA)) + n1*log((n1*nA)/(s1 + n1*nA))))/(v*(s1 + n1*nA)) + (chi*kT*nA*((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1)*(n - n1))/(V - v*(s1 + n1*nA)) + (chi*kT*n1*nA*s1)/(v*(s1 + n1*nA)^2)) + (2473475804651973*gamma*v)/(70368744177664*(V - v*(s1 + n1*nA))^(1/3)))/kT) - beta + k2*n1*nA + (alpha*n1*nA)/(v*(s1 + n1*nA))^(2/3) - k1*nA*r*v*(s1 + n1*nA) - (alpha*n1*nA*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT))/(v*(s1 + n1*nA))^(2/3);
dX(3) = V*k1*r - k2*n;
dX(4) = c1 - c2*r;

A(1, 1) = (alpha*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT))/(v*(s1 + n1*nA))^(2/3) - alpha/(v*(s1 + n1*nA))^(2/3) - k2 + k1*nA*r*v + (2*alpha*n1*nA*v)/(3*(v*(s1 + n1*nA))^(5/3)) - (2*alpha*n1*nA*v*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT))/(3*(v*(s1 + n1*nA))^(5/3)) + (alpha*n1*nA*v*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT)*((kT*((nA*((nA*v)/(V - v*(s1 + n1*nA)) - (nA^2*v^2*(n - n1))/(V - v*(s1 + n1*nA))^2))/((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1) + (2*chi*nA^2*v)/(V - v*(s1 + n1*nA)) - (((nA*v)/(V - v*(s1 + n1*nA)) - (nA^2*v^2*(n - n1))/(V - v*(s1 + n1*nA))^2)*(V - v*(s1 + n1*nA)))/(nA*v*(n - n1)) - (2*chi*nA^3*v^2*(n - n1))/(V - v*(s1 + n1*nA))^2))/(nA*v) - (kT*(nA^2/(s1 + n1*nA) - (2*chi*nA^2)/(s1 + n1*nA) + (2*chi*n1*nA^3)/(s1 + n1*nA)^2 + ((s1 + n1*nA)*(nA/(s1 + n1*nA) - (n1*nA^2)/(s1 + n1*nA)^2))/(n1*nA)))/(nA*v)))/(kT*(v*(s1 + n1*nA))^(2/3));
A(1, 2) = k1*r*v + (2*alpha*n1*v)/(3*(v*(s1 + n1*nA))^(5/3)) - (2*alpha*n1*v*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT))/(3*(v*(s1 + n1*nA))^(5/3)) - (alpha*n1*nA*v*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT)*((kT*((nA*(s1 + n1*nA)*(s1/(s1 + n1*nA)^2 - 1/(s1 + n1*nA)))/s1 - 1/(s1 + n1*nA) + (2*chi*n1*nA^2)/(s1 + n1*nA)^2))/(nA*v) + (kT*((nA^2*v^2*(n - n1))/(((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1)*(V - v*(s1 + n1*nA))^2) - v/(V - v*(s1 + n1*nA)) + (2*chi*nA^2*v^2*(n - n1))/(V - v*(s1 + n1*nA))^2))/(nA*v)))/(kT*(v*(s1 + n1*nA))^(2/3));
A(1, 3) = -(alpha*n1*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT)*((nA^2*v)/(((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1)*(V - v*(s1 + n1*nA))) - 1/(n - n1) + (2*chi*nA^2*v)/(V - v*(s1 + n1*nA))))/(v*(s1 + n1*nA))^(2/3);
A(1, 4) = k1*v*(s1 + n1*nA);
A(2, 1) = k2*nA + (alpha*nA)/(v*(s1 + n1*nA))^(2/3) - k1*nA^2*r*v - (alpha*nA*exp((nA*v*((kT*(log(-(nA*v*(n - n1))/(s1*v - V + n1*nA*v)) - nA - nA*log((s1*v - V + n*nA*v)/(s1*v - V + n1*nA*v)) + chi*nA + (2*chi*nA^2*v*(n - n1))/(s1*v - V + n1*nA*v) + 1))/(nA*v) + (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT))/(v*(s1 + n1*nA))^(2/3) - (2*alpha*n1*nA^2*v)/(3*(v*(s1 + n1*nA))^(5/3)) - (beta*exp((v*((kT*(n - n1)*(log(-(nA*v*(n - n1))/(s1*v - V + n1*nA*v)) - nA - nA*log((s1*v - V + n*nA*v)/(s1*v - V + n1*nA*v)) + chi*nA + (2*chi*nA^2*v*(n - n1))/(s1*v - V + n1*nA*v) + 1))/(s1*v - V + n1*nA*v) - (kT*n1*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(v*(s1 + n1*nA))) - v*((kT*(s1*log(s1/(s1 + n1*nA)) + n1*log((n1*nA)/(s1 + n1*nA))))/(v*(s1 + n1*nA)) - (kT*((log((s1*v - V + n*nA*v)/(s1*v - V + n1*nA*v))*(s1*v - V + n*nA*v))/(s1*v - V + n1*nA*v) - (v*log(-(nA*v*(n - n1))/(s1*v - V + n1*nA*v))*(n - n1))/(s1*v - V + n1*nA*v)))/v + (chi*kT*nA*(n - n1)*(s1*v - V + n*nA*v))/(s1*v - V + n1*nA*v)^2 + (chi*kT*n1*nA*s1)/(v*(s1 + n1*nA)^2)) + (2473475804651973*gamma*v)/(70368744177664*(V - s1*v - n1*nA*v)^(1/3)))/kT)*(70368744177664*V^3*kT*s1^2 + 70368744177664*V^3*kT*n1^2*nA^3 - 140737488355328*V^2*kT*s1^3*v + 70368744177664*V*kT*s1^4*v^2 + 70368744177664*V^3*kT*n1*nA^2*s1 + 824491934883991*gamma*n1^3*nA^4*v^2*(V - s1*v - n1*nA*v)^(5/3) + 70368744177664*kT*n*nA*s1^4*v^3 + 70368744177664*V^2*kT*n1^3*nA^3*v - 211106232532992*V^2*kT*n1^3*nA^4*v - 70368744177664*V*kT*n1^4*nA^4*v^2 + 140737488355328*V*kT*n1^4*nA^5*v^2 + 70368744177664*kT*n*n1^4*nA^5*v^3 - 70368744177664*kT*n*n1^4*nA^6*v^3 - 70368744177664*kT*n*nA^2*s1^4*v^3 + 70368744177664*V^3*kT*n1*nA*s1 + 824491934883991*gamma*nA*s1^3*v^2*(V - s1*v - n1*nA*v)^(5/3) + 70368744177664*V*kT*n*nA^2*s1^3*v^2 - 211106232532992*V^2*kT*n1*nA^2*s1^2*v + 140737488355328*V*kT*n1*nA^2*s1^3*v^2 - 422212465065984*V^2*kT*n1^2*nA^3*s1*v - 140737488355328*V*kT*n1^3*nA^3*s1*v^2 + 422212465065984*V*kT*n1^3*nA^4*s1*v^2 - 140737488355328*chi*kT*n*n1^4*nA^6*v^3 + 140737488355328*chi*kT*n*nA^2*s1^4*v^3 + 281474976710656*kT*n*n1*nA^2*s1^3*v^3 - 281474976710656*kT*n*n1*nA^3*s1^3*v^3 + 281474976710656*kT*n*n1^3*nA^4*s1*v^3 - 281474976710656*kT*n*n1^3*nA^5*s1*v^3 + 422212465065984*V*kT*n1^2*nA^3*s1^2*v^2 + 140737488355328*chi*kT*n^2*n1^3*nA^6*v^3 + 140737488355328*chi*kT*n^2*nA^3*s1^3*v^3 + 422212465065984*kT*n*n1^2*nA^3*s1^2*v^3 - 422212465065984*kT*n*n1^2*nA^4*s1^2*v^3 - 140737488355328*V^3*chi*kT*n1*nA^2*s1 - 70368744177664*V*kT*n*nA*s1^3*v^2 - 211106232532992*V^2*kT*n1*nA*s1^2*v + 140737488355328*V*kT*n1*nA*s1^3*v^2 + 2473475804651973*gamma*n1*nA^2*s1^2*v^2*(V - s1*v - n1*nA*v)^(5/3) + 2473475804651973*gamma*n1^2*nA^3*s1*v^2*(V - s1*v - n1*nA*v)^(5/3) + 140737488355328*V*chi*kT*n1^4*nA^5*v^2 - 70368744177664*V*kT*n*n1^3*nA^4*v^2 + 70368744177664*V*kT*n*n1^3*nA^5*v^2 - 140737488355328*V*chi*kT*n*n1^3*nA^5*v^2 - 140737488355328*V*chi*kT*n*nA^2*s1^3*v^2 + 422212465065984*V^2*chi*kT*n1*nA^2*s1^2*v - 281474976710656*V*chi*kT*n1*nA^2*s1^3*v^2 + 422212465065984*V^2*chi*kT*n1^2*nA^3*s1*v - 211106232532992*V*kT*n*n1*nA^2*s1^2*v^2 + 211106232532992*V*kT*n*n1*nA^3*s1^2*v^2 - 211106232532992*V*kT*n*n1^2*nA^3*s1*v^2 + 211106232532992*V*kT*n*n1^2*nA^4*s1*v^2 + 281474976710656*chi*kT*n*n1*nA^3*s1^3*v^3 - 281474976710656*chi*kT*n*n1^3*nA^5*s1*v^3 - 422212465065984*V*chi*kT*n1^2*nA^3*s1^2*v^2 + 422212465065984*chi*kT*n^2*n1*nA^4*s1^2*v^3 + 422212465065984*chi*kT*n^2*n1^2*nA^5*s1*v^3 - 422212465065984*V*chi*kT*n*n1*nA^3*s1^2*v^2 - 422212465065984*V*chi*kT*n*n1^2*nA^4*s1*v^2))/(70368744177664*kT*(s1 + n1*nA)^3*(s1*v - V + n1*nA*v)^3) + (2*alpha*n1*nA^2*v*exp((nA*v*((kT*(log(-(nA*v*(n - n1))/(s1*v - V + n1*nA*v)) - nA - nA*log((s1*v - V + n*nA*v)/(s1*v - V + n1*nA*v)) + chi*nA + (2*chi*nA^2*v*(n - n1))/(s1*v - V + n1*nA*v) + 1))/(nA*v) + (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT))/(3*(v*(s1 + n1*nA))^(5/3)) + (alpha*n1*nA^2*v*exp((nA*v*((kT*(log(-(nA*v*(n - n1))/(s1*v - V + n1*nA*v)) - nA - nA*log((s1*v - V + n*nA*v)/(s1*v - V + n1*nA*v)) + chi*nA + (2*chi*nA^2*v*(n - n1))/(s1*v - V + n1*nA*v) + 1))/(nA*v) + (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT)*((kT*(s1^2 + n1^2*nA^3 + n1*nA^2*s1 + n1*nA*s1 - 2*chi*n1*nA^2*s1))/(n1*nA*v*(s1 + n1*nA)^2) + (kT*(s1^2*v^2 + V^2 - 2*V*s1*v + n1^2*nA^3*v^2 - V*n*nA*v - V*n1*nA*v + 2*chi*n^2*nA^3*v^2 + V*n*nA^2*v - V*n1*nA^2*v + n*nA*s1*v^2 + n1*nA*s1*v^2 + n*n1*nA^2*v^2 - n*n1*nA^3*v^2 - n*nA^2*s1*v^2 + n1*nA^2*s1*v^2 - 2*chi*n*n1*nA^3*v^2 + 2*chi*n*nA^2*s1*v^2 - 2*chi*n1*nA^2*s1*v^2 - 2*V*chi*n*nA^2*v + 2*V*chi*n1*nA^2*v))/(nA*v*(n - n1)*(s1*v - V + n1*nA*v)^2)))/(kT*(v*(s1 + n1*nA))^(2/3));
A(2, 2) = (beta*exp((v*((kT*(n - n1)*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(V - v*(s1 + n1*nA)) - (kT*n1*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(v*(s1 + n1*nA))) - v*((kT*(log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA)))*((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1) - (v*log((nA*v*(n - n1))/(V - v*(s1 + n1*nA)))*(n - n1))/(V - v*(s1 + n1*nA))))/v + (kT*(s1*log(s1/(s1 + n1*nA)) + n1*log((n1*nA)/(s1 + n1*nA))))/(v*(s1 + n1*nA)) + (chi*kT*nA*((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1)*(n - n1))/(V - v*(s1 + n1*nA)) + (chi*kT*n1*nA*s1)/(v*(s1 + n1*nA)^2)) + (2473475804651973*gamma*v)/(70368744177664*(V - v*(s1 + n1*nA))^(1/3)))/kT)*(v*((kT*((v^2*(n - n1))/(V - v*(s1 + n1*nA))^2 + (v^2*log((nA*v*(n - n1))/(V - v*(s1 + n1*nA)))*(n - n1))/(V - v*(s1 + n1*nA))^2 - (nA*v^2*(n - n1))/(V - v*(s1 + n1*nA))^2 - (nA*v^2*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA)))*(n - n1))/(V - v*(s1 + n1*nA))^2))/v + (kT*(n1/(s1 + n1*nA) - log(s1/(s1 + n1*nA)) + (s1 + n1*nA)*(s1/(s1 + n1*nA)^2 - 1/(s1 + n1*nA))))/(v*(s1 + n1*nA)) + (kT*(s1*log(s1/(s1 + n1*nA)) + n1*log((n1*nA)/(s1 + n1*nA))))/(v*(s1 + n1*nA)^2) - (chi*kT*n1*nA)/(v*(s1 + n1*nA)^2) - (chi*kT*nA^2*v^2*(n - n1)^2)/(V - v*(s1 + n1*nA))^3 + (2*chi*kT*n1*nA*s1)/(v*(s1 + n1*nA)^3) - (chi*kT*nA*v*((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1)*(n - n1))/(V - v*(s1 + n1*nA))^2) + v*((kT*(n - n1)*((nA^2*v^2*(n - n1))/(((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1)*(V - v*(s1 + n1*nA))^2) - v/(V - v*(s1 + n1*nA)) + (2*chi*nA^2*v^2*(n - n1))/(V - v*(s1 + n1*nA))^2))/(V - v*(s1 + n1*nA)) + (kT*n1*((nA*(s1 + n1*nA)*(s1/(s1 + n1*nA)^2 - 1/(s1 + n1*nA)))/s1 - 1/(s1 + n1*nA) + (2*chi*n1*nA^2)/(s1 + n1*nA)^2))/(v*(s1 + n1*nA)) + (kT*v*(n - n1)*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(V - v*(s1 + n1*nA))^2 + (kT*n1*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(v*(s1 + n1*nA)^2)) + (824491934883991*gamma*v^2)/(70368744177664*(V - v*(s1 + n1*nA))^(4/3))))/kT - k1*nA*r*v - (2*alpha*n1*nA*v)/(3*(v*(s1 + n1*nA))^(5/3)) + (2*alpha*n1*nA*v*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT))/(3*(v*(s1 + n1*nA))^(5/3)) + (alpha*n1*nA^2*v*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT)*((kT*((nA*(s1 + n1*nA)*(s1/(s1 + n1*nA)^2 - 1/(s1 + n1*nA)))/s1 - 1/(s1 + n1*nA) + (2*chi*n1*nA^2)/(s1 + n1*nA)^2))/(nA*v) + (kT*((nA^2*v^2*(n - n1))/(((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1)*(V - v*(s1 + n1*nA))^2) - v/(V - v*(s1 + n1*nA)) + (2*chi*nA^2*v^2*(n - n1))/(V - v*(s1 + n1*nA))^2))/(nA*v)))/(kT*(v*(s1 + n1*nA))^(2/3));
A(2, 3) = (beta*exp((v*((kT*(n - n1)*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(V - v*(s1 + n1*nA)) - (kT*n1*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(v*(s1 + n1*nA))) - v*((kT*(log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA)))*((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1) - (v*log((nA*v*(n - n1))/(V - v*(s1 + n1*nA)))*(n - n1))/(V - v*(s1 + n1*nA))))/v + (kT*(s1*log(s1/(s1 + n1*nA)) + n1*log((n1*nA)/(s1 + n1*nA))))/(v*(s1 + n1*nA)) + (chi*kT*nA*((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1)*(n - n1))/(V - v*(s1 + n1*nA)) + (chi*kT*n1*nA*s1)/(v*(s1 + n1*nA)^2)) + (2473475804651973*gamma*v)/(70368744177664*(V - v*(s1 + n1*nA))^(1/3)))/kT)*(v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(V - v*(s1 + n1*nA)) + (kT*(n - n1)*((nA^2*v)/(((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1)*(V - v*(s1 + n1*nA))) - 1/(n - n1) + (2*chi*nA^2*v)/(V - v*(s1 + n1*nA))))/(V - v*(s1 + n1*nA))) - v*((chi*kT*nA*((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(V - v*(s1 + n1*nA)) - (kT*(v/(V - v*(s1 + n1*nA)) + (v*log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))))/(V - v*(s1 + n1*nA)) - (nA*v)/(V - v*(s1 + n1*nA)) - (nA*v*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))))/(V - v*(s1 + n1*nA))))/v + (chi*kT*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA))^2)))/kT + (alpha*n1*nA*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT)*((nA^2*v)/(((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1)*(V - v*(s1 + n1*nA))) - 1/(n - n1) + (2*chi*nA^2*v)/(V - v*(s1 + n1*nA))))/(v*(s1 + n1*nA))^(2/3);
A(2, 4) = -k1*nA*v*(s1 + n1*nA);
A(3, 1) = 0;
A(3, 2) = 0;
A(3, 3) = -k2;
A(3, 4) = V*k1;
A(4, 1) = 0;
A(4, 2) = 0;
A(4, 3) = 0;
A(4, 4) = -c2;

BBT(1, 1) = k2*n1 + (alpha*n1)/(v*(s1 + n1*nA))^(2/3) + (alpha*n1*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT))/(v*(s1 + n1*nA))^(2/3) + k1*r*v*(s1 + n1*nA);
BBT(1, 2) = - k2*n1*nA - (alpha*n1*nA)/(v*(s1 + n1*nA))^(2/3) - k1*nA*r*v*(s1 + n1*nA) - (alpha*n1*nA*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT))/(v*(s1 + n1*nA))^(2/3);
BBT(1, 3) = k2*n1 + k1*r*v*(s1 + n1*nA);
BBT(1, 4) = 0;
BBT(2, 1) = - k2*n1*nA - (alpha*n1*nA)/(v*(s1 + n1*nA))^(2/3) - k1*nA*r*v*(s1 + n1*nA) - (alpha*n1*nA*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT))/(v*(s1 + n1*nA))^(2/3);
BBT(2, 2) = beta + beta*exp((v*((kT*(n - n1)*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(V - v*(s1 + n1*nA)) - (kT*n1*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(v*(s1 + n1*nA))) - v*((kT*(log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA)))*((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1) - (v*log((nA*v*(n - n1))/(V - v*(s1 + n1*nA)))*(n - n1))/(V - v*(s1 + n1*nA))))/v + (kT*(s1*log(s1/(s1 + n1*nA)) + n1*log((n1*nA)/(s1 + n1*nA))))/(v*(s1 + n1*nA)) + (chi*kT*nA*((nA*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1)*(n - n1))/(V - v*(s1 + n1*nA)) + (chi*kT*n1*nA*s1)/(v*(s1 + n1*nA)^2)) + (2473475804651973*gamma*v)/(70368744177664*(V - v*(s1 + n1*nA))^(1/3)))/kT) + k2*n1*nA^2 + (alpha*n1*nA^2)/(v*(s1 + n1*nA))^(2/3) + k1*nA^2*r*v*(s1 + n1*nA) + (alpha*n1*nA^2*exp(-(nA*v*((kT*(nA - log((nA*v*(n - n1))/(V - v*(s1 + n1*nA))) - chi*nA + nA*log(1 - (nA*v*(n - n1))/(V - v*(s1 + n1*nA))) + (2*chi*nA^2*v*(n - n1))/(V - v*(s1 + n1*nA)) - 1))/(nA*v) - (kT*(nA - log((n1*nA)/(s1 + n1*nA)) - chi*nA + nA*log(s1/(s1 + n1*nA)) + (2*chi*n1*nA^2)/(s1 + n1*nA) - 1))/(nA*v)))/kT))/(v*(s1 + n1*nA))^(2/3);
BBT(2, 3) = - k2*n1*nA - k1*nA*r*v*(s1 + n1*nA);
BBT(2, 4) = 0;
BBT(3, 1) = k2*n1 + k1*r*v*(s1 + n1*nA);
BBT(3, 2) = - k2*n1*nA - k1*nA*r*v*(s1 + n1*nA);
BBT(3, 3) = k2*n + V*k1*r;
BBT(3, 4) = 0;
BBT(4, 1) = 0;
BBT(4, 2) = 0;
BBT(4, 3) = 0;
BBT(4, 4) = c1 + c2*r;

dC = A*C + C*A' + BBT;
dN = [dX(:); reshape(dC, 16, 1)];
