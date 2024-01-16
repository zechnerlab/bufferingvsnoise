function NMat = CalculateConcentrationDependence_Ternary_NonEQ_Correlation(phiTotGridA, phiTotGridB, grid, chiAA, chiAB, chiBB, kT, v, nA, V, gamma, alpha, beta, k2A, k2B, c1, c2, plotTraces)

if nargin<18
    plotTraces = 0;
end

numStates = 6;
mur = c1/c2; % mean RNA abundance

NMat = zeros(length(phiTotGridA), length(phiTotGridB), numStates + numStates^2);

for j=1:length(phiTotGridA)
    for u=1:length(phiTotGridB)
        phiTotA = phiTotGridA(j);
        phiTotB = phiTotGridB(u);
        a = V*phiTotA/(nA*v);
        b = V*phiTotB/(nA*v);
        s = (V - (a+b)*nA*v)/v;
        k1A = a*k2A/(mur*V);
        k1B = b*k2B/(mur*V);
        
        % Set initial conditions
        N0 = zeros(numStates + numStates^2, 1);
        N0(1) = a-1; %initialize almost everything in dilute phase and a small droplet with volume fraction 0.5.
        N0(2) = b-1;
        N0(3) = s-1;
        N0(4) = a;
        N0(5) = b;
        N0(6) = mur;
        
        % Solve LNA differential equations generated above
        options = {};
        [gridF, N] = ode15s(@ActiveTernaryMixture_Correlation, grid, N0, options, kT, v, chiAA, chiAB, chiBB, V, nA, alpha, beta, gamma, k1A, k1B, k2A, k2B, c1, c2);
        
        % Plotting
        if (plotTraces == 1)
            figure(1);
            subplot(2,3,1);
            plot(gridF, N(:, 1), gridF, N(:, 1) - sqrt(N(:, 7)), gridF, N(:, 1) + sqrt(N(:, 7)));
            
            subplot(2,3,2);
            plot(gridF, N(:, 1)*nA./((N(:, 1) + N(:, 2))*nA + N(:, 3)));
            
            subplot(2,3,3);
            plot(gridF, (N(:, 4) - N(:, 1))*nA*v./(V - v*((N(:, 1) + N(:, 2))*nA + N(:, 3))));
            
            subplot(2,3,4);
            plot(gridF, V - v*((N(:, 1) + N(:, 2))*nA + N(:, 3)));
            
            subplot(2,3,5);
            plot(gridF, N(:, 4), gridF, N(:, 4) + sqrt(N(:, 28)), gridF, N(:, 3) - sqrt(N(:, 28)));
            
            subplot(2,3,6);
            plot(gridF, N(:, 5));
        end
        
        % Store end points of the moment matrix
        NMat(j, u, :) = N(end, :);
        
        % Some printing
        fprintf('Finished iteration %d (%d); %d (%d).\n', j, u, length(phiTotGridA), length(phiTotGridB));
        
    end
end