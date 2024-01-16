%% Calculate statistics of phi1 and phiTot.
% Moments of volume fractions are calculated using a Taylor expansion.
% Comparison to Gaussian integration is provided below but commented out.

function dropletStats = CalculateStatisticsBinary_NonEQ_SDL(NMat, nA, v, V, numStates, phi2)

dropletStats.N1_Mu = NMat(:, 1);
dropletStats.N1_Var = NMat(:, 4);
dropletStats.N1NTot_Cov = NMat(:, 5);

dropletStats.NTot = NMat(:, 2);
dropletStats.NTot_Var = NMat(:, 8);

dropletStats.phi_mu = dropletStats.NTot*nA*v/V;
dropletStats.phi_Var = dropletStats.NTot_Var*(nA*v/V)^2;


dropletStats.phi1_mu = dropletStats.N1_Mu*nA*v./V;%(V - nA*v*(dropletStats.NTot - dropletStats.N1_Mu)/phi2);
dropletStats.phi1_Var = dropletStats.N1_Var*(nA*v/V)^2;

%% Store Mu and Sigma matrix for each condition and compare moments obtained 
% by Taylor expansion with Gaussian integral (latter one commented out)

for k=1:size(NMat, 1)
    dropletStats.Mu{k} = NMat(k, 1:numStates)';
    dropletStats.Sigma{k} = reshape(NMat(k, numStates+1:end), numStates, numStates);
end



end