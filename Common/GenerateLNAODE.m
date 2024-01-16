function GenerateLNAODE(S, h, params, stateNames, functionName)
numStates = size(S, 2);
if nargin<4
    X = sym('X', [numStates, 1]);
else
    X = sym(stateNames);
end
if nargin<5
   functionName = 'LNAMomentsODE';
end
fileName = [functionName '.m'];

h = sym(h(:));
A = simplify(jacobian(transpose(S)*h, X));
BBT = simplify(transpose(S)*diag(h)*S);
F = simplify(transpose(S)*h);


fid = fopen(fileName, 'w');

fprintf(fid, 'function dN = %s(t, N', functionName);

for k=1:length(params)
   fprintf(fid, ', %s', params{k});
end
fprintf(fid, ')\n');

for k=1:numStates
    fprintf(fid, '%s=N(%d);\n', char(X(k)), k);
end

fprintf(fid, 'C = reshape(N(%d:end), %d, %d);\n\n', numStates+1, numStates, numStates);

for k=1:numStates
fprintf(fid, 'dX(%d) = %s;\n', k, char(F(k, :)));
end

fprintf(fid, '\n');

for k=1:numStates
   for j=1:numStates
    fprintf(fid, 'A(%d, %d) = %s;\n', k, j, char(A(k, j)));
   end
end

fprintf(fid, '\n');

for k=1:numStates
   for j=1:numStates
    fprintf(fid, 'BBT(%d, %d) = %s;\n', k, j, char(BBT(k, j)));
   end
end
   
fprintf(fid, '\n');
fprintf(fid, 'dC = A*C + C*A'' + BBT;');
fprintf(fid, '\n');
fprintf(fid, 'dN = [dX(:); reshape(dC, %d, 1)];\n', numStates*numStates);

fclose(fid);