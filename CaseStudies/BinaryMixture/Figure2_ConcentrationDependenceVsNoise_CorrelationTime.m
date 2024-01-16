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
load results/results_Figure2.mat;

%% Simulate the generated system using the linear noise approximation

%This value of phiTot corresponds to the second ellipse in the concentration-dependency plot
phiTot = results.dropletStats.phi_mu(max(results.plotIdx));

c2Grid = logspace(-7, 0, 25);
c2Grid = sort([c2Grid, c2]);
markExampleIdx = find(c2Grid == c2); %needed to plot the point from the other figure panel showing the ellipses.
muR = 10;


%% Calculate phiPlus under the dilute approximation (zero osmotic pressure)
fun = @(phiPlus) OsmoticPressureBinary_FH_Approx(phiPlus, chi, nA);

options = optimset('Display', 'off', 'TolX', 1e-18);
phi2_Approx = fzero(fun, 0.6, options);
f2 = FreeEnergyDensityBinary_FH(phi2_Approx, chi, kT, v, nA);
nStar = (exp(-1 + nA - nA*chi + (nA*v*f2)/(kT*phi2_Approx))*V)/(nA*v);
phiStar = nStar*nA*v / V;

dropletStats_SDL = {};
dropletStats = {};

for k=1:length(c2Grid)
    c2 = c2Grid(k);
    c1 = muR*c2;
    
    % Sweep phi total to obtain concentration dependence and fluctuations at
    % equilibrium    
    NMat = CalculateConcentrationDependence_Binary_SDL(phiTot, grid, nStar, v, nA, V, kD, k2, c1, c2, plotTraces);
    
    % Calculate moments for volume fractions
    dropletStats_SDL{k} = CalculateStatisticsBinary_NonEQ_SDL(NMat, nA, v, V, 3, phi2_Approx);
    
    
    % Sweep phi total to obtain concentration dependence and fluctuations at
    % equilibrium
    %phiTotGrid = linspace(0.0003, 0.01, 50);
    [NMat] = CalculateConcentrationDependence_Binary_NonEQ(phiTot, grid, chi, kT, v, nA, V, gamma, alpha, beta, k2, c1, c2, plotTraces);
    
    % Calculate moments for volume fractions
    dropletStats{k} = CalculateStatisticsBinary_NonEQ(NMat, nA, v, V, 4);
    
    
end

%% Calculate noise reduction for different c2

plotSamples = 1;
M = 1000000;
numSamplesShown = 20; %number of samples
NoiseReduction = zeros(size(c2Grid));
NoiseReduction_SDL = zeros(size(c2Grid));

% define colormap
cmap = [linspace(1, 0.6, 256)', linspace(1, 0.8, 256)', linspace(1, 0.9, 256)'];

% For each selected condition, read out Mu and Sigma and draw M samples.
% These samples are used to (a) create a contour plot of the joint
% distribution over phi and phi1 and (b) to draw individual samples.
for k = 1:length(c2Grid)
    
    Mu = dropletStats{k}.Mu{1};
    Sigma = dropletStats{k}.Sigma{1};
    
    X = Mu + sqrtm(Sigma)*randn(numStates, M);
    phi1_s = X(1, :)*nA./(X(1, :)*nA + X(2, :));
    phiTot_s = X(3, :)*nA*v/V;
    
    NoiseReduction(k) = ((dropletStats{k}.phi1_Var(1)/dropletStats{k}.phi1_mu(1)^2) / (dropletStats{k}.phi_Var(1)/dropletStats{k}.phi_mu(1)^2))^(-1/2);
    NoiseReduction_SDL(k) = ((dropletStats_SDL{k}.phi1_Var(1)/dropletStats_SDL{k}.phi1_mu(1)^2) / (dropletStats_SDL{k}.phi_Var(1)/dropletStats_SDL{k}.phi_mu(1)^2))^(-1/2);
    fprintf('NoiseTot/Noise1 = %f (c2=%f)\n', NoiseReduction(k), c2Grid(k));
    
end


bufferingStrength = 1 + phiStar*kD./(k2*phiTot);

loglog((1./c2Grid)/60, NoiseReduction, 'Color', cols(2, :), 'LineWidth', 2); hold on;
loglog((1./c2Grid(markExampleIdx))/60, NoiseReduction(markExampleIdx), 'o', 'Color', cols(2, :), 'MarkerFaceColor', cols(2, :));
loglog((1./c2Grid)/60, NoiseReduction_SDL, '--', 'Color', cols(2, :), 'LineWidth', 2); 
loglog((1./c2Grid)/60, repmat(bufferingStrength, 1, length(c2Grid)), '--', 'Color', cols(1, :), 'LineWidth', 2); 
xlabel('Time in min');
ylabel('Noise Reduction');
ylim([1e-1, 2e1]);
xlim([1e-1, 1e5]);
xticks([0.1, 10, 1000, 100000]);
box off;
currAx = gca;
currAx.FontSize = 10;
currAx.XColor = 'k';
currAx.YColor = 'k';
