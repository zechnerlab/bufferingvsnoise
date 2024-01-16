function [NMat, k1Vec] = CalculateConcentrationDependence_Binary_NonEQ(phiTotGrid, grid, chi, kT, v, nA, V, gamma, alpha, beta, k2, c1, c2, plotTraces)

if nargin<14
   plotTraces = 0; 
end

numStates = 4;
mur = c1/c2; % mean RNA abundance

NMat = zeros(length(phiTotGrid), numStates + numStates^2);

for u=1:length(phiTotGrid)
    phiTot = phiTotGrid(u);
    n = V*phiTot/(nA*v);
    s = (V - n*nA*v)/v;
    k1 = n*k2/(mur*V);
    k1Vec(u) = k1;
    
    % Set initial conditions
    N0 = zeros(numStates + numStates^2, 1);
    N0(1) = n-1; %initialize almost everything in dilute phase and a small droplet with volume fraction 0.5.
    N0(2) = s-1;
    N0(3) = n;
    N0(4) = mur;
    
    % Solve LNA differential equations generated above
    options = {};
    [gridF, N] = ode15s(@ActiveBinaryMixture, grid, N0, options, kT, v, chi, V, nA, alpha, beta, gamma, k1, k2, c1, c2);
    
    % Plotting
    if (plotTraces == 1)
        figure(1);
        subplot(2,3,1);
        plot(gridF, N(:, 1), gridF, N(:, 1) - sqrt(N(:, 5)), gridF, N(:, 1) + sqrt(N(:, 5)));
        
        subplot(2,3,2);
        plot(gridF, N(:, 1)*nA./(N(:, 1)*nA + N(:, 2)));
        
        subplot(2,3,3);
        plot(gridF, (N(:, 3) - N(:, 1))*nA*v./(V - v*(N(:, 1)*nA + N(:, 2))));
        
        subplot(2,3,4);
        plot(gridF, (V - v*(N(:, 1)*nA + N(:, 2))));
        
        subplot(2,3,5);
        plot(gridF, N(:, 3), gridF, N(:, 3) + sqrt(N(:, 20)), gridF, N(:, 3) - sqrt(N(:, 20)));
        
        subplot(2,3,6);
        plot(gridF, N(:, 4));
    end
    
    % Store end points of the moment matrix 
    NMat(u, :) = N(end, :);

    % Some printing
    fprintf('Finished iteration %d of %d (phiTot=%f).\n', u, length(phiTotGrid), phiTot);
    
end