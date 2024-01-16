clear;
close all;

addpath('../../../Common/Graphics/')

rng(1000);

cols = GetDefaultColors();

%% get steady state data
data = readtable('Klosin2020_2NTDDX4_Fig2c.xlsx');
ConcNucYFPnM   = table2array(data(:,{'total_mean_conc'}));   % intensity of YFP averaged across the whole nucleus
ConcPlsmYFPnM  = table2array(data(:,{'plasm_mean_conc'}));    % intensity of YFP outside of the droplets
AreaNucleus  = table2array(data(:,{'nuclear_area'}));    % area total nucleus
AreaDroplets = table2array(data(:,{'Z_drops_Area'}));    % area total droplet phase


[sortA, sortIdx] = sort(AreaDroplets);

% choose cells that are well in the phase coexistence regime (i.e., c>12uM)
selIdx = ConcNucYFPnM>12;
X = ConcNucYFPnM(selIdx);
Y = ConcPlsmYFPnM(selIdx);

% choosen number of intervals were 34 and 68
K = 34;
qGrid = linspace(0, 1, K+1);
qX = quantile(X, qGrid);

for k=1:length(qX)-1
   selIdx = and(X> qX(k), X<=qX(k+1));
   XSel = X(selIdx);
   YSel = Y(selIdx);
   
   numSamples = length(XSel);
   n_bootstrap = 10000;
   rdIdx = randi(numSamples, numSamples, n_bootstrap);
   XSel_resample = XSel(rdIdx);
   YSel_resample = YSel(rdIdx);
   
   meanTot_bootstrap = mean(XSel_resample);
   meanDilute_bootstrap = mean(YSel_resample);
   meanTot(k) = mean(meanTot_bootstrap);
   meanDilute(k) = mean(meanDilute_bootstrap);
   
   sigmaMat = cov(meanTot_bootstrap, meanDilute_bootstrap);
   varMeanTot(k) = sigmaMat(1, 1);
   varMeanDilute(k) = sigmaMat(2, 2);
   covMeanTotMeanDilute(k) = sigmaMat(1, 2);

end

errorbar(meanTot, meanDilute, 2.5*sqrt(varMeanDilute), 2.5*sqrt(varMeanDilute), 2.5*sqrt(varMeanTot), 2.5*sqrt(varMeanTot), '.');
hold on;


%% Infer phiStar and kD/k2 using MCMC
M = 100000;
propSigma = 0.1;
chain = zeros(2, M);
chain(:, 1) = log([5, 5]);

for k=2:M

    proposal = chain(:, k-1);
    
    for i=1:length(proposal)
        proposal(i) = normrnd(chain(i, k-1), propSigma);
    end
    
    [LNew, phi1Pred] = evaluateLogLikelihood(meanTot, meanDilute, exp(proposal(1)), exp(proposal(2)), varMeanTot, varMeanDilute, covMeanTotMeanDilute);
    
    if (k>2)
       a = min(1, exp(LNew - LOld)); 
    else
       a = 1;
    end
    
    if (rand < a)
       %fprintf('Accepted\n');
       chain(:, k) = proposal;
       LOld = LNew;
       currPred = phi1Pred;
    else
        %fprintf('Rejected\n');
        chain(:, k) = chain(:, k-1);
       
    end
    
    if (mod(k, 2000)==2)
        
        subplot(1,2,1);
        plot(exp(chain(:, 1:k)')); 
        xlabel('Iterations');
        ylabel('Parameter values');
        
        subplot(1,2,2);
        plot(meanTot, meanDilute, 'o', meanTot, currPred, 'x');
        xlabel('\langle c \rangle (\muM)');
        ylabel('\langle c_1 \rangle (\muM)');
        
        drawnow;
        
        fprintf('Finishd iteration %d (%d).\n', k, M);
    end
    
end

% Extract estimates and uncertainties
chainLinearDom = exp(chain);
burnIn = 10000;
phiStar_MCMC = mean(chainLinearDom(1, burnIn:end));
phiStar_MCMC_Std = std(chainLinearDom(1, burnIn:end));
kD_k2_Ratio_MCMC = mean(chainLinearDom(2, burnIn:end));
kD_k2_Ratio_MCMC_Std = std(chainLinearDom(2, burnIn:end));

fprintf('Phi_sat = %f +/- %f\n', phiStar_MCMC, phiStar_MCMC_Std);
fprintf('kD/k2 = %f +/- %f\n', kD_k2_Ratio_MCMC, kD_k2_Ratio_MCMC_Std); 

% Select all cells with total concentration larger than 10nM to estimate
% phiStar and kD/k2 using a simple linear regression for comparison.
X_linReg = ConcNucYFPnM(ConcNucYFPnM>10);
Y_linReg  = ConcPlsmYFPnM(ConcNucYFPnM>10);

modelFit = fitlm(X_linReg, Y_linReg);

a = modelFit.Coefficients.Estimate(1); %offset
b = modelFit.Coefficients.Estimate(2); %slope

phiStar_linReg = a / (1-b);
kD_k2_Ratio_linReg = a/phiStar_linReg / b;

totConcGrid = linspace(0, 120, 100);
concPlsmPred = gApprox(totConcGrid, phiStar_MCMC, kD_k2_Ratio_MCMC);

figure;
plot(ConcNucYFPnM, ConcPlsmYFPnM, '.', 'Color', [0.85, 0.5, 0.4], 'MarkerSize', 8);hold on;
plot([0, 120], [phiStar_MCMC, phiStar_MCMC], 'k--', [phiStar_MCMC, phiStar_MCMC], [0, 30], 'k--', 'LineWidth', 2);
plot(totConcGrid, concPlsmPred, '--', 'Color', cols(1, :), 'LineWidth', 2);
plot([0, phiStar_MCMC], [0, phiStar_MCMC], '--', 'Color', cols(1, :), 'LineWidth', 2);
plot(meanTot, meanDilute, '.', 'Color', cols(1, :), 'MarkerSize', 18);
box off;
xlabel('c (\muM)');
ylabel('c_1 (\muM)');
ylim([0, 30]);
xlim([0, 120]);
currAx = gca;
currAx.FontSize = 12;
currAx.XColor = 'k';
currAx.YColor = 'k';

data = readtable('Klosin2020_2NTDDX4_dense_dilute_FigureS22.csv');
concDrop_smallset  = table2array(data(:,{'dense'}));    % area total nucleus
concPlsm_smallset = table2array(data(:,{'dilute'}));

figure;
plot(concPlsm_smallset, concDrop_smallset ./ concPlsm_smallset, '.', 'Color', [0.85, 0.5, 0.4], 'MarkerSize', 8); hold on;

meanDropletConc = mean(concDrop_smallset);

concPlsmTheory = linspace(min(concPlsm_smallset), max(concPlsm_smallset), 100);
rhoTheory = meanDropletConc ./ concPlsmTheory;


plot(concPlsmTheory, rhoTheory, '--', 'Color', cols(1, :), 'LineWidth', 2)
xlabel('c_1 (\muM)');
ylabel('c_2 / c_1');
box off;
ylim([50, 200]);
yticks([50, 100, 150, 200]);
xlim([5, 15]);
xticks([5, 10, 15]);
currAx = gca;
currAx.FontSize = 10;
currAx.XColor = 'k';
currAx.YColor = 'k';

fileName = sprintf('results/Inference_K=%d.mat', K);
save(fileName);

