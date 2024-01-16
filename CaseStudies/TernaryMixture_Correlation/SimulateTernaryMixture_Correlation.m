%% This script analyzes fluctuations in binary phase separating systems with
% active molecular turnover. Solutions are calculated with the linear noise
% approximation and compared to the analytical solution obtained in the
% small droplet limit.
% Author: Christoph Zechner (zechner@mpi-cbg.de)

clear;
close all;

addpath('../../Common');
addpath('../../TernaryMixture');
addpath('../../TernaryMixture_Correlation');

cols = GetDefaultColors();

% Set to one to output some traces of the dynamics (for visualization / debugging).
plotTraces = 0;

%% Load parameter configuration. Check 'CreateConfig.m' for details.
regionB = 2;
filename = sprintf('data/Config_regionB=%d.mat', regionB);
load(filename);

%% Calculate concentration-dependence at thermodynamic equilibrium 
% using ternary Flory-Huggins model.
[phiA1EQ, phiB1EQ, phiA2EQ, phiB2EQ, V2EQ] = CalculateConcentrationDependence_TernaryFH(phiTotGridAEQ, phiTotB, chiAA, chiAB, chiBB, kT, v, nA, V, gamma);


%% Calculate concentration-dependence and fluctuations away from equilibrium
[NMat] = CalculateConcentrationDependence_Ternary_NonEQ_Correlation(phiTotGridA, ...
    phiTotB, grid, chiAA, chiAB, chiBB, kT, v, nA, V, gamma, alpha, beta, k2A, k2B, c1, c2, plotTraces);

% Calculate moments for volume fractions
dropletStats = CalculateStatisticsTernary_NonEQ_Correlation(NMat, nA, v, V, numStates);

%% Save results.
fileNameOut = sprintf('data/Results_Figure2_Config_regionB=%d.mat', regionB);
save(fileNameOut);
