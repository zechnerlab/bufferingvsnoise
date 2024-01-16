%% This script analyzes fluctuations in ternary phase separating systems with
% active molecular turnover. 
% Author: Christoph Zechner (zechner@mpi-cbg.de)

clear;
close all;

% initialize seed for random number generation
rng(1000);

addpath('../../Common');
addpath('../../TernaryMixture');

cols = GetDefaultColors();

%% Load results created by "SimulateTernaryMixture.m"

regionBIdx = [2];

for u=1:length(regionBIdx)
    
    regionB = regionBIdx(u);
    fileName = sprintf('data/Results_Figure2_Config_regionB=%d.mat', regionB);
    load(fileName);
    
    %% Draw statistics of individual condition(s)
    plotIdxVec = ceil([0.2, 0.5]*length(phiTotGridA));%fliplr(ceil([0.3]*length(phiTotGridA)));
    markExampleIdx = max(plotIdxVec); %index of the condition that is shown (here the upper of the 2 ellipses).
    
    % decide whether individual random samples of phi1 and phiTot should be
    % plotted
    plotSamples = 1;
    M = 1000000;
    numSamplesShown = 20; %number of samples
    NoiseReductionA = zeros(size(plotIdxVec));
    NoiseReductionB = zeros(size(plotIdxVec));
    
    % define colormap
    cmap = [linspace(1, 0.9, 256)', linspace(1, 0.3, 256)', linspace(1, 0.3, 256)'];
    
    % For each selected condition, read out Mu and Sigma and draw M samples.
    % These samples are used to (a) create a contour plot of the joint
    % distribution over phi and phi1 and (b) to draw individual samples.
    
    for k = length(phiTotGridA):-1:1 %plot reverse because otherwise the background of the larger ellipse will overlap with the lower one
        
        Mu = dropletStats.Mu{k, 1};
        Sigma = dropletStats.Sigma{k, 1};
        
        X = Mu + sqrtm(Sigma)*randn(numStates, M);
        phiA1_s = X(1, :)*nA./((X(1, :) + X(2, :))*nA + X(3, :));
        phiA2_s = (X(4, :) - X(1, :))*nA*v./(V - ((X(1, :) + X(2, :))*nA + X(3, :))*v);
        phiTotA_s = X(4, :)*nA*v/V;
        phiB1_s = X(2, :)*nA./((X(1, :) + X(2, :))*nA + X(3, :));
        phiTotB_s = X(5, :)*nA*v/V;
        phiB2_s = (X(5, :) - X(2, :))*nA*v./(V - ((X(1, :) + X(2, :))*nA + X(3, :))*v);
        
        NoiseReductionA1(k) = (var(phiA1_s) / mean(phiA1_s)^2 / (var(phiTotA_s) / mean(phiTotA_s)^2))^(-1/2);
        NoiseReductionB1(k) = (var(phiB1_s) / mean(phiB1_s)^2 / (var(phiTotB_s) / mean(phiTotB_s)^2))^(-1/2);
        NoiseReductionA2(k) = (var(phiA2_s) / mean(phiA2_s)^2 / (var(phiTotA_s) / mean(phiTotA_s)^2))^(-1/2);
        NoiseReductionB2(k) = (var(phiB2_s) / mean(phiB2_s)^2 / (var(phiTotB_s) / mean(phiTotB_s)^2))^(-1/2);
        
        
        
        idx = find(k==plotIdxVec);
        
        if (~isempty(idx))
            plotIdx = plotIdxVec(idx);
            
            figure(10+u);
            [Z, edgesPhiATot, edgesPhiA1] = histcounts2(phiTotA_s, phiA1_s, 'Normalization', 'probability');
            [~, c] = contourf(edgesPhiATot(1:end-1), edgesPhiA1(1:end-1), Z', 'LevelStep', 50e-5);
            c.LineWidth = 0.5;
            c.LineColor = [0.7, 0.7, 0.7];
            colormap(cmap);
            hold on;
            
            figure(20+u);
            subplot(1, length(plotIdxVec), idx);
            [Z, edgesPhiATot, edgesPhiA1] = histcounts2(phiTotA_s/mean(phiTotA_s), phiA1_s/mean(phiA1_s), 'Normalization', 'probability');
            [~, c] = contourf(edgesPhiATot(1:end-1), edgesPhiA1(1:end-1), Z', 'LevelStep', 50e-5);
            c.LineWidth = 0.5;
            c.LineColor = [0.7, 0.7, 0.7];
            colormap(cmap);
            hold on;
            xlabel('\phi / \langle \phi \rangle');
            ylabel('\phi_1 / \langle \phi_1 \rangle');
            
            
            if (plotSamples == 1)
                figure(10+u);
                plot(phiTotA_s(1:numSamplesShown), phiA1_s(1:numSamplesShown), '.', 'Color', cols(2, :)); hold on;
                
                figure(20+u);
                subplot(1, length(plotIdxVec), idx);
                plot(phiTotA_s(1:numSamplesShown)/mean(phiTotA_s), phiA1_s(1:numSamplesShown)/mean(phiA1_s), '.', 'Color', cols(2, :)); hold on;
                xlim([0.75, 1.25]);
                ylim([0.75, 1.25]);
                box off;
                
            end
            
            fprintf('Noise1/NoiseTot (component A) = %f (phi=%f)\n', NoiseReductionA1(plotIdx), phiTotGridA(plotIdx));
            
        end 
    end
    
    figure(10+u);
    plot([phiTotGridA], [dropletStats.phiA1_mu], '-','LineWidth', 2, 'Color', cols(1, :)); hold on;
    plot(phiTotGridAEQ, phiA1EQ, '-', 'LineWidth', 1.5, 'Color', 'k');
    xlabel('Total volume fraction \phi');
    ylabel('Dilute phase volume fraction \phi^I');
    ylim([0, 4e-4]);
    xlim([0, 0.03]);
    box off;
    currAx = gca;
    currAx.FontSize = 12;
    currAx.XColor = 'k';
    currAx.YColor = 'k';
    
    dLogPhi = diff(log(phiTotGridA));
    dLogPhi1 = diff(log(dropletStats.phiA1_mu))';
    
    figure;
    semilogy(phiTotGridA, NoiseReductionA1, 'Color', cols(2, :), 'LineWidth', 2); hold on;
    semilogy(phiTotGridA(1:end-1), abs(dLogPhi./dLogPhi1), 'Color', cols(1, :), 'LineWidth', 2); 
    xlabel('Total volume fraction \langle \phi \rangle');
    ylabel('Noise reduction');
    box off;
    legend('\Gamma', '\Lambda');
    currAx = gca;
    currAx.FontSize = 10;
    currAx.XColor = 'k';
    currAx.YColor = 'k';
    ylim([0.1, 1000]);
    yticks([0.1, 10, 1000]);
    xlim([0, 0.03]);
    xticks([0, 0.01, 0.02, 0.03]);
    
end

