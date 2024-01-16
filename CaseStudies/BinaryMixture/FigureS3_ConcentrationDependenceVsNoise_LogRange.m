%% This script analyzes fluctuations in binary phase separating systems with
% active molecular turnover. Solutions are calculated with the linear noise
% approximation and compared to the analytical solution obtained in the
% small droplet limit.
% Author: Christoph Zechner (zechner@mpi-cbg.de)

clear;
close all;

addpath('../../Common');
addpath('../../BinaryMixture');
addpath('../../BinaryMixture_ZOPA');

cols = GetDefaultColors();

% Set flags
plotTraces = 0;

% load results / settings from main analysis.
load results/results_FigureS3.mat;

%% Calculate concentration-dependence at thermodynamic equilibrium 
% using binary Flory-Huggins model.
phiTotMax = 0.02;
phiTotGridEQ = linspace(0.00001, phiTotMax, 200);
[phi1_mu_eq, phi1_Var_eq] = CalculateConcentrationDependence_BinaryFH(phiTotGridEQ, chi, kT, v, nA, V, gamma);

%% Calculate phiPlus under the dilute approximation (zero osmotic pressure)
fun = @(phiPlus) OsmoticPressureBinary_FH_Approx(phiPlus, chi, nA);

options = optimset('Display', 'off', 'TolX', 1e-18);
phi2_Approx = fzero(fun, 0.6, options);
f2 = FreeEnergyDensityBinary_FH(phi2_Approx, chi, kT, v, nA);
nStar = (exp(-1 + nA - nA*chi + (nA*v*f2)/(kT*phi2_Approx))*V)/(nA*v);
phiStar = nStar*nA*v / V;

% Sweep phi total to obtain concentration dependence and fluctuations at
% equilibrium
phiTotExample = results.dropletStats.phi_mu(max(results.plotIdx));
phiTotGrid = logspace(-3, log10(phiTotMax), 100);
phiTotGrid = sort([phiTotGrid, phiTotExample]);
markExampleIdx = find(phiTotGrid==phiTotExample);

NMat = CalculateConcentrationDependence_Binary_SDL(phiTotGrid, grid, nStar, v, nA, V, kD, k2, c1, c2, plotTraces);

% Calculate moments for volume fractions
dropletStats_SDL = CalculateStatisticsBinary_NonEQ_SDL(NMat, nA, v, V, 3, phi2_Approx);


% Sweep phi total to obtain concentration dependence and fluctuations at
% equilibrium
%phiTotGrid = linspace(0.0003, 0.01, 50);
[NMat] = CalculateConcentrationDependence_Binary_NonEQ(phiTotGrid, grid, chi, kT, v, nA, V, gamma, alpha, beta, k2, c1, c2, plotTraces);

% Calculate moments for volume fractions
dropletStats = CalculateStatisticsBinary_NonEQ(NMat, nA, v, V, 4);


%% Draw statistics of individual condition(s)
NoiseReduction = zeros(size(phiTotGrid));
NoiseReduction_SDL = zeros(size(phiTotGrid));
M = 1000000;

% For each selected condition, read out Mu and Sigma and draw M samples.
% These samples are used to (a) create a contour plot of the joint
% distribution over phi and phi1 and (b) to draw individual samples.
subplot(1,2,1);
for k = 1:length(phiTotGrid)
    
    
    Mu = dropletStats.Mu{k};%NMat(plotIdx, 1:numStates)';
    Sigma = dropletStats.Sigma{k};%reshape(NMat(plotIdx, numStates+1:end), numStates, numStates);
    
    X = Mu + sqrtm(Sigma)*randn(numStates, M);
    phi1_s = X(1, :)*nA./(X(1, :)*nA + X(2, :));
    phiTot_s = X(3, :)*nA*v/V;
    
    NoiseReduction(k) = ((dropletStats.phi1_Var(k)/dropletStats.phi1_mu(k)^2) / (dropletStats.phi_Var(k)/dropletStats.phi_mu(k)^2))^(-1/2);
    NoiseReduction_SDL(k) = ((dropletStats_SDL.phi1_Var(k)/dropletStats_SDL.phi1_mu(k)^2) / (dropletStats_SDL.phi_Var(k)/dropletStats_SDL.phi_mu(k)^2))^(-1/2);
    fprintf('NoiseTot/Noise1 = %f (phi=%f)\n', NoiseReduction(k), phiTotGrid(k));
    
end

%% Estimate noise reduction from concentration dependence g and compare to true noise reduction
dLogPhi = diff(log(phiTotGrid));
dLogPhiPlus = diff(log(dropletStats.phi1_mu))';
figure;
semilogy(phiTotGrid(1:end-1), dLogPhi./dLogPhiPlus, 'LineWidth', 2, 'Color', cols(1, :)); hold on;
semilogy(phiTotGrid, 1 + phiStar*kD./(k2*phiTotGrid), '--', 'LineWidth', 2, 'Color', cols(1, :));
semilogy(phiTotGrid, NoiseReduction, 'LineWidth',2, 'Color', cols(2, :));
semilogy(phiTotGrid, NoiseReduction_SDL, '--', 'LineWidth', 2, 'Color', cols(2, :));
xlabel('\langle \phi \rangle');
ylabel('Buffering strength / noise reduction');
box off;
legend('\Lambda', '\Lambda (dil. approx.)', '\Gamma', '\Gamma (dil. approx.)');
currAx = gca;
currAx.FontSize = 12;
currAx.XColor = 'k';
currAx.YColor = 'k';
xticks([0, 0.005, 0.01, 0.015, 0.02]);


