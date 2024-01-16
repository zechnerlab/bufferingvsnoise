%% This script analyzes fluctuations in binary phase separating systems with
% active molecular turnover. Solutions are calculated with the linear noise
% approximation and compared to the analytical solution obtained in the
% small droplet limit.
% Author: Christoph Zechner (zechner@mpi-cbg.de)

clear;
close all;

% initialize seed for random number generation
rng(1000);

addpath('../../Common');
addpath('../../BinaryMixture');
addpath('../../BinaryMixture_DA');

cols = GetDefaultColors();

% Set flags
plotTraces = 0;


%% Simulate the generated system using the linear noise approximation
numStates = 4;
T = 1000000*60;
grid = linspace(0, T, 100);

kT = 1;
% Molecular volume v
v = 1/20;
% Nondimensional molecular volume of the solute nA. If nA=2, for instance,
% the solute is twice the size of the solvent.
nA = 20;
% System Volume V
V = 10e6;
% Flory-Huggins interaction parameter chi
chi = 1.3;
% Surface tension
gamma = 1e-5;%0*1e-1;% 0*4e-2 ;

% Fast point: chi=6, V = 500000, nA=1.2, alpha = 20, beta = 1000, phiTot
% between 0.001-0.3 or so.`

% Diffusion parameter (subsuming factor 6)
alpha = 1000;
% Permeability (sort of)
beta = 500000;
% kD
kD = alpha / (V^(2/3));

% Birth and death rate constants for protein
k1 = 0.1; % this one is irrelevant since it will be sweeped later.
k2 = 0.000005; % death rate

% Birth and death rate constants for RNA
c1 = 0.003; % birth (default: 0.001)
c2 = 0.0005; % death (default: 0.0003);


%% Calculate concentration-dependence at thermodynamic equilibrium 
% using binary Flory-Huggins model.
phiTotMax = 0.005;
phiTotGridEQ = linspace(0.000001, phiTotMax, 100);
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
phiTotGrid = linspace(phiStar*1.3, phiTotMax, 20);
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
plotIdxVec = fliplr(ceil([0.3, 0.7]*length(phiTotGrid)));
lineStyles = {'-', '-'};

% decide whether individual random samples of phi1 and phiTot should be
% plotted
plotSamples = 1;
M = 1000000;
numSamplesShown = 20; %number of samples
NoiseReduction = zeros(size(plotIdxVec));


% define colormap
cmap = [linspace(1, 0.9, 256)', linspace(1, 0.3, 256)', linspace(1, 0.3, 256)'];

% For each selected condition, read out Mu and Sigma and draw M samples.
% These samples are used to (a) create a contour plot of the joint
% distribution over phi and phi1 and (b) to draw individual samples.
%subplot(1,2,1);

for k = 1:length(plotIdxVec)
    plotIdx = plotIdxVec(k);
    
    Mu = dropletStats.Mu{plotIdx};
    Sigma = dropletStats.Sigma{plotIdx};
    
    X = Mu + sqrtm(Sigma)*randn(numStates, M);
    phi1_s = X(1, :)*nA./(X(1, :)*nA + X(2, :));
    phiTot_s = X(3, :)*nA*v/V;
    
    figure(10);
    [Z, edgesPhiTot, edgesPhi1] = histcounts2(phiTot_s, phi1_s, 'Normalization', 'probability');
    [~, c] = contourf(edgesPhiTot(1:end-1), edgesPhi1(1:end-1), Z', 'LevelStep', 50e-5);
    c.LineWidth = 0.5;
    c.LineColor = [0.7, 0.7, 0.7];
    c.LineStyle = lineStyles{k};
    colormap(cmap); 
    hold on;
    currAx = gca;
    currAx.FontSize = 10;
    currAx.XColor = 'k';
    currAx.YColor = 'k';
    
    figure(10+k);
    [Z, edgesPhiTot, edgesPhi1] = histcounts2(phiTot_s/mean(phiTot_s), phi1_s/mean(phi1_s), 'Normalization', 'probability');
    [~, c] = contourf(edgesPhiTot(1:end-1), edgesPhi1(1:end-1), Z', 'LevelStep', 50e-5);
    c.LineWidth = 0.5;
    c.LineColor = [0.7, 0.7, 0.7];
    c.LineStyle = lineStyles{k};
    colormap(cmap); 
    xlim([0.75, 1.25 ]);
    ylim([0.75, 1.25 ]);
    box off;
    hold on;
    currAx = gca;
    currAx.FontSize = 10;
    currAx.XColor = 'k';
    currAx.YColor = 'k';
    xlabel('\phi');
    ylabel('\phi_1');
    
    Noise1(k) = sqrt(dropletStats.phi1_Var(plotIdx)/dropletStats.phi1_mu(plotIdx)^2);
    NoiseTot(k) = sqrt(dropletStats.phi_Var(plotIdx)/dropletStats.phi_mu(plotIdx)^2);
    NoiseReduction(k) = ((dropletStats.phi1_Var(plotIdx)/dropletStats.phi1_mu(plotIdx)^2) / (dropletStats.phi_Var(plotIdx)/dropletStats.phi_mu(plotIdx)^2))^(-1/2);
    fprintf('NoiseTot=%f, Noise1=%f, NoiseTot/Noise1 = %f (phi=%f)\n', NoiseTot(k), Noise1(k), NoiseReduction(k), phiTotGrid(plotIdx));
    
    if (plotSamples == 1)
        figure(10);
        plot(phiTot_s(1:numSamplesShown), phi1_s(1:numSamplesShown), '.', 'Color', cols(2, :));
    
        figure(10+k);
        plot(phiTot_s(1:numSamplesShown)/mean(phiTot_s), phi1_s(1:numSamplesShown)/mean(phi1_s), '.', 'Color', cols(2, :));
    end
end

figure(10);
plot([phiTotGrid], [dropletStats.phi1_mu], 'LineWidth', 2, 'Color', cols(1, :)); hold on;
plot([phiTotGrid], [dropletStats_SDL.phi1_mu], '--', 'LineWidth', 2, 'Color', cols(1, :));
plot(phiTotGridEQ, phi1_mu_eq, '-', 'LineWidth', 1.5, 'Color', 'k');
xlabel('Total volume fraction \phi');
ylabel('Dilute phase volume fraction \phi^I');
xlim([0, max(phiTotGrid)]);
ylim([0, 0.6e-4]);
box off;
currAx = gca;
currAx.FontSize = 12;
currAx.XColor = 'k';
currAx.YColor = 'k';
xticks([0, 0.001, 0.002, 0.003, 0.004, 0.005]);


results.NoiseReduction = NoiseReduction;
results.dropletStats = dropletStats;
results.plotIdx = plotIdxVec;
save results/results_FigureS3.mat;
