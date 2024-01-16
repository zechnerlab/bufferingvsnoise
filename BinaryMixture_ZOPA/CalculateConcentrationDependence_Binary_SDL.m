function NMat = CalculateConcentrationDependence_Binary_SDL(phiTotGrid, grid, nStar, v, nA, V, kD, k2, c1, c2, plotTraces)

if nargin<11
   plotTraces = 0; 
end

mur = c1/c2; % mean RNA abundance

numStates = 3;

NMat = zeros(length(phiTotGrid), numStates + numStates^2);

for u=1:length(phiTotGrid)
    phiTot = phiTotGrid(u);
    n = V*phiTot/(nA*v);
    k1 = n*k2/(mur*V);
    
    % Set initial conditions
    N0 = zeros(numStates + numStates^2, 1);
    N0(1) = n-1; %initialize almost everything in dilute phase and a small droplet with volume fraction 0.5.
    N0(2) = n;
    N0(3) = mur;
    
    % Solve LNA differential equations generated above
    options = {};
    [gridF, N] = ode15s(@ActiveBinaryMixture_SmallDropletLimit, grid, N0, options, V, kD, nStar, k1, k2, c1, c2);
    
    % Plotting
    if (plotTraces == 1)
        figure(1);
        subplot(1,3,1);
        plot(gridF, N(:, 1), gridF, N(:, 1) - sqrt(N(:,4)), gridF, N(:, 1) + sqrt(N(:, 4)));
        
        subplot(1,3,2);
        plot(gridF, N(:, 1)*nA*v./V);
 
        subplot(1,3,3);
        plot(gridF, N(:, 2), gridF, N(:, 2) + sqrt(N(:, 8)), gridF, N(:, 2) - sqrt(N(:, 8)));
    end
    
    % Store end points of the moment matrix 
    NMat(u, :) = N(end, :);

    % Some printing
    fprintf('Finished iteration %d of %d (phiTot=%f).\n', u, length(phiTotGrid), phiTot);
    
end
