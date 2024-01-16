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


%% Simulate the generated system using the linear noise approximation
numStates = 4;
T = 1000000*60;
grid = linspace(0, T, 100);

kT = 1;
% Molecular volume v
v = 1/100;
% Nondimensional molecular volume of the solute nA. If nA=2, for instance,
% the solute is twice the size of the solvent.
nA = 20;
% System Volume V
V = 10e6;
% Flory-Huggins interaction parameter chi
chi = 1.2;
% Surface tension
gamma = 1e-5;%0*1e-1;% 0*4e-2 ;

% Fast point: chi=6, V = 500000, nA=1.2, alpha = 20, beta = 1000, phiTot
% between 0.001-0.3 or so.`

% Timescale ratio vector kD/k2
timeScaleRatio = [200, 1000];

% Permeability (sort of)
beta = 500000;

% Birth and death rate constants for protein
k1 = 0.1; % this one is irrelevant since it will be sweeped later.
k2 = 0.000005; % death rate

% Birth and death rate constants for RNA
c1 = 0.001; % birth (default: 0.001)
c2 = 0.0001; % death (default: 0.0003);


%% Calculate concentration-dependence at thermodynamic equilibrium
% using binary Flory-Huggins model.
phiTotMax = 0.02;
phiTotGridEQ = linspace(0.00001, phiTotMax, 100);
[phi1_mu_eq, phi1_Var_eq, phi2_mu_eq, ~, V2_eq] = CalculateConcentrationDependence_BinaryFH(phiTotGridEQ, chi, kT, v, nA, V, gamma);


%% Calculate phiPlus under the dilute approximation (zero osmotic pressure)
fun = @(phiPlus) OsmoticPressureBinary_FH_Approx(phiPlus, chi, nA);

options = optimset('Display', 'off', 'TolX', 1e-18);
phi2_Approx = fzero(fun, 0.6, options);
f2 = FreeEnergyDensityBinary_FH(phi2_Approx, chi, kT, v, nA);
nStar = (exp(-1 + nA - nA*chi + (nA*v*f2)/(kT*phi2_Approx))*V)/(nA*v);
phiStar = nStar*nA*v / V;


figure(10);
plot(phiTotGridEQ, phi1_mu_eq, '-', 'LineWidth', 2, 'Color', 'k'); hold on;

xlabel('Total volume fraction \phi');
ylabel('Dilute phase volume fraction \phi^I');
box off;

figure(11);
twoPhaseIdx = phiTotGridEQ>phiStar*1.3; %Define threshold to focus on regime of phase coexistence.
plot(phiTotGridEQ(twoPhaseIdx), phi2_mu_eq(twoPhaseIdx)./phi1_mu_eq(twoPhaseIdx), '-', 'LineWidth', 2, 'Color', 'k'); hold on;
xlabel('Total volume fraction \phi');
ylabel('Partition coefficient');
box off;

for k=1:length(timeScaleRatio)
    
    kD = timeScaleRatio(k)*k2;
    alpha = kD * V^(2/3);
    
    % Sweep phi total to obtain concentration dependence and fluctuations at
    % equilibrium
    phiTotGrid = linspace(phiStar*1.3, phiTotMax, 20); %Focus on two phase regime (same threshold as above)
    NMat = CalculateConcentrationDependence_Binary_SDL(phiTotGrid, grid, nStar, v, nA, V, kD, k2, c1, c2, plotTraces);
    
    % Calculate moments for volume fractions
    dropletStats_SDL = CalculateStatisticsBinary_NonEQ_SDL(NMat, nA, v, V, 3, phi2_Approx);
    
    % Sweep phi total to obtain concentration dependence and fluctuations at
    % equilibrium
    %phiTotGrid = linspace(0.0003, 0.01, 50);
    [NMat] = CalculateConcentrationDependence_Binary_NonEQ(phiTotGrid, grid, chi, kT, v, nA, V, gamma, alpha, beta, k2, c1, c2, plotTraces);
    
    % Calculate moments for volume fractions
    dropletStats = CalculateStatisticsBinary_NonEQ(NMat, nA, v, V, 4);
    
    figure(10);
    plot([phiTotGrid], [dropletStats.phi1_mu], 'LineWidth', 2, 'Color', cols(1, :)); hold on;
    plot([phiTotGrid], [dropletStats_SDL.phi1_mu], '--', 'LineWidth', 2, 'Color', cols(1, :));
    
    figure(11);
    plot([phiTotGrid], [dropletStats.phi2_mu./dropletStats.phi1_mu], 'LineWidth', 2, 'Color', cols(1, :)); hold on;
    plot([phiTotGrid], [phi2_Approx./dropletStats_SDL.phi1_mu], '--', 'LineWidth', 2, 'Color', cols(1, :));
    drawnow;
    
end
xlim([0, max(phiTotGrid)]);
%ylim([0, 3e-4]);

