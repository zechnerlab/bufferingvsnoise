%% Calculate statistics of phi1 and phiTot.
% Moments of volume fractions are calculated using a Taylor expansion.
% Comparison to Gaussian integration is provided below but commented out.

function dropletStats = CalculateStatisticsBinary_NonEQ(NMat, nA, v, V, numStates)

dropletStats.N1_Mu = NMat(:, 1);
dropletStats.N1_Var = NMat(:, 5);
dropletStats.S1_Mu = NMat(:, 2);
dropletStats.S1_Var = NMat(:, 10);
dropletStats.N1S1_Cov = NMat(:, 6);
dropletStats.N1NTot_Cov = NMat(:, 7);
dropletStats.S1NTot_Cov = NMat(:, 11);

dropletStats.NTot = NMat(:, 3);
dropletStats.NTot_Var = NMat(:, 15);

dropletStats.phi_mu = dropletStats.NTot*nA*v/V;
dropletStats.phi_Var = dropletStats.NTot_Var*(nA*v/V)^2;


dropletStats.phi1_mu = dropletStats.N1_Mu*nA./(dropletStats.N1_Mu*nA + dropletStats.S1_Mu);
dropletStats.phi2_mu = (dropletStats.NTot - dropletStats.N1_Mu)*nA*v./(V - v*(dropletStats.N1_Mu*nA + dropletStats.S1_Mu));
dropletStats.V2 = (V - v*(dropletStats.N1_Mu*nA + dropletStats.S1_Mu));


dropletStats.phi1_Var = (nA^2*(dropletStats.S1_Var.*dropletStats.N1_Mu.^2 + ...
    dropletStats.S1_Mu.*(-2*dropletStats.N1S1_Cov.*dropletStats.N1_Mu + ...
    dropletStats.N1_Var.*dropletStats.S1_Mu)))./(nA*dropletStats.N1_Mu + dropletStats.S1_Mu).^4;
dropletStats.V2_Var = v^2*(-2*dropletStats.N1S1_Cov*nA + nA^2*dropletStats.N1_Var + dropletStats.S1_Var);

dropletStats.phiTotphi1_Cov = (nA^2*v*(-dropletStats.S1NTot_Cov.*dropletStats.N1_Mu + dropletStats.N1NTot_Cov.*dropletStats.S1_Mu ...
   + dropletStats.NTot.*dropletStats.N1_Mu.*(nA*dropletStats.N1_Mu + ...
   dropletStats.S1_Mu)))./(V*(nA*dropletStats.N1_Mu + dropletStats.S1_Mu).^2) - dropletStats.phi1_mu.*dropletStats.phi_mu;

dropletStats.phiTotphi1_Cov = (nA^2*v*(-dropletStats.S1NTot_Cov.*dropletStats.N1_Mu + dropletStats.N1NTot_Cov.*dropletStats.S1_Mu))...
    ./(V*(nA*dropletStats.N1_Mu + dropletStats.S1_Mu).^2);


%% Store Mu and Sigma matrix for each condition and compare moments obtained 
% by Taylor expansion with Gaussian integral (latter one commented out)

for k=1:size(NMat, 1)
    dropletStats.Mu{k} = NMat(k, 1:numStates)';
    dropletStats.Sigma{k} = reshape(NMat(k, numStates+1:end), numStates, numStates);
%     
%     X = Mu + sqrtm(Sigma)*randn(numStates, 1000);
%     phi1_s = X(1, :)*nA./(X(1, :)*nA + X(2, :));
%     phiTot_s = X(3, :)*nA/V;
%     
%     %subplot(2,2,1); hold on;
%     %plot(phiTot_s, phi1_s, '.');
%     
%     Mu0 = Mu(1:2);
%     Sigma0 = Sigma(1:2, 1:2); 
%     
%     funMean = @(n, s) GaussianIntegral_VolumeFraction(n, s, nA, 1, Mu0, Sigma0);
%     funSecondMoment = @(n, s) GaussianIntegral_VolumeFraction(n, s, nA, 2, Mu0, Sigma0);
%     
%     sf = 5;
%     phi1_mu2(k) = integral2(funMean, Mu0(1)-sf*sqrt(Sigma0(1, 1)), ...
%         Mu0(1)+sf*sqrt(Sigma0(1, 1)), Mu0(2)-sf*sqrt(Sigma0(2, 2)), Mu0(2)+sf*sqrt(Sigma0(2, 2)));
%     phi1_var2(k) = integral2(funSecondMoment, Mu0(1)-sf*sqrt(Sigma0(1, 1)), ...
%         Mu0(1)+sf*sqrt(Sigma0(1, 1)), Mu0(2)-sf*sqrt(Sigma0(2, 2)), Mu0(2)+sf*sqrt(Sigma0(2, 2)))...
%         -phi1_mu2(k)^2;
end



end