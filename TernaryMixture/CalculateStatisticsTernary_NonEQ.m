%% Calculate statistics of phi1 and phiTot.
% Calculates mean vectors and covariance matrices and stores them in a
% structure dropletStats

function dropletStats = CalculateStatisticsTernary_NonEQ(NMat, nA, v, V, numStates)


dropletStats.A1_Mu = NMat(:, :, 1);
dropletStats.A1_Var = NMat(:, :, 8);
dropletStats.B1_Mu = NMat(:, :, 2);
dropletStats.B1_Var = NMat(:, :, 16);
dropletStats.S1_Mu = NMat(:, :, 3);
dropletStats.S1_Var = NMat(:, :, 24);
dropletStats.ATot_Mu = NMat(:, :, 4);
dropletStats.BTot_Mu = NMat(:, :, 5);

dropletStats.phiA1_mu = nA*dropletStats.A1_Mu ./ ((dropletStats.A1_Mu + dropletStats.B1_Mu)*nA + dropletStats.S1_Mu);
dropletStats.phiATot = dropletStats.ATot_Mu*nA*v/V;
dropletStats.phiB1_mu = nA*dropletStats.B1_Mu ./ ((dropletStats.A1_Mu + dropletStats.B1_Mu)*nA + dropletStats.S1_Mu);
dropletStats.phiBTot = dropletStats.BTot_Mu*nA*v/V;

dropletStats.phiA2_mu = nA*v*(dropletStats.ATot_Mu - dropletStats.A1_Mu)./ (V - v*((dropletStats.A1_Mu + dropletStats.B1_Mu)*nA + dropletStats.S1_Mu));
dropletStats.phiB2_mu = nA*v*(dropletStats.BTot_Mu - dropletStats.B1_Mu)./ (V - v*((dropletStats.A1_Mu + dropletStats.B1_Mu)*nA + dropletStats.S1_Mu));

%% Store Mu and Sigma matrix for each condition and compare moments obtained 
% by Taylor expansion with Gaussian integral (latter one commented out)

numRows = size(NMat, 1);
numCols = size(NMat, 2);

for j=1:numRows
    for k=1:numCols
        mVec = NMat(j, k, 1:numStates);
        dropletStats.Mu{j, k} = mVec(:);
        
        sigmaVec = NMat(j, k, numStates+1:end);
        dropletStats.Sigma{j, k} = reshape(sigmaVec(:), numStates, numStates);
    end
end


end