
%% Create config files for simulations
%regionB sets the value of phiBTot (can be 1, 2, 3)
clear;
regionB = 2;

numStates = 6;
T = 10000*60;
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
chiAA = 0;
chiAB = -4;
chiBB = 0;

%chiAA = 0.8;
%chiAB = -2;
%chiBB = 0;

% Surface tension
gamma = 1e-5;

% Diffusion parameter (subsuming factor 6)
alpha = 500;
% Permeability (sort of)
beta = 500000;
% kD
kD = alpha / (V^(2/3));

% Birth and death rate constants for protein
k1A = 0.1; % this one is irrelevant since it will be sweeped later.
k1B = 0.1;
k2A = 0.00002; % death rate A
k2B = 0.00002; % death rate B

% Birth and death rate constants for RNA
c1 = 0.001; % birth (default: 0.001)
c2 = 0.0001; % death (default: 0.0003);

% Set values for sweep in phiTotA / phiTotB
phiTotMax = 0.03;
% Different grids are used for equilibrium / non-equilibrium due to
% different runtimes / numerics
phiTotGridAEQ = linspace(0.00001, phiTotMax, 300);

% Grid for simulating the whole phase diagram
phiTotGridAPD = linspace(0.00001, 0.4, 25);
phiTotGridBPD = linspace(0.00001, 0.4, 25);
% Threshold needed to filter out conditions where there is no stable
% droplet (i.e., to get only the points corresponding to phase
% coexistence). Since the optimization can give small but non-zero droplet
% volumes, a threshold (in percent of the total volume V) needs to be
% provided.
PD_VolumeThreshold = 0.001;

numPointsNonEQ = 50;

if (regionB==1)
   phiTotGridA = linspace(0.01, phiTotMax, 10);
   phiTotB = 0.01;
elseif (regionB==2)
   phiTotGridA = linspace(0.0007, phiTotMax, numPointsNonEQ);
   phiTotB = 0.045;
elseif (regionB==3)
   phiTotGridA = linspace(0.0007, 0.04, numPointsNonEQ);
   phiTotB = 0.06; 
elseif (regionB==4)
   phiTotGridA = linspace(0.002, 0.05, 10);
   phiTotB = 0.06; 
end

fileName = sprintf('data/Config_regionB=%d.mat', regionB);
save(fileName);

