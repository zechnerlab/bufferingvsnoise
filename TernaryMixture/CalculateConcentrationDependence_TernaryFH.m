function [phiA1Opt, phiB1Opt, phiA2Opt, phiB2Opt, V2Opt, phiTotAMat, phiTotBMat] = CalculateConcentrationDependence_TernaryFH(phiTotGridA, phiTotGridB, chiAA, chiAB, chiBB, kT, v, nA, V, gamma)

% Initialize vectors for the calculated solutions
phiA1Opt = zeros(length(phiTotGridA), length(phiTotGridB));
phiB1Opt = phiA1Opt;
V2Opt = phiA1Opt;
phiA2Opt = phiA1Opt;
phiB2Opt = phiB1Opt;
phiTotAMat = phiA1Opt;
phiTotBMat = phiB1Opt;

%Intitialize optimizer
options = optimoptions('fmincon', 'Algorithm', 'interior-point');
options.Display = 'None';

fprintf('\nCalculating equilibrium concentration-dependence of binary FH model...\n');

for k=1:length(phiTotGridA)
    for j=1:length(phiTotGridB)
        % Set total volume fraction in the system and calculate how many total
        % solutes and solvents this corresponds to
        phiTotA = phiTotGridA(k);
        phiTotB = phiTotGridB(j);
        phiTotAMat(k, j) = phiTotA;
        phiTotBMat(k, j) = phiTotB;
        
        a = V*phiTotA/(nA*v);
        b = V*phiTotB/(nA*v);
        s = (V - (a + b)*nA*v)/v;
        
        % Set function handle to FH Free energy
        fun = @(state) FreeEnergyTernary_FH(state(1), state(2), state(3), chiAA, chiAB, chiBB, kT, v, nA, V, a, b, gamma);
        
        % TODO: try multistart optimization to find multiple optima.
        % Find minimum of the Free Energy
        
        % initialize a smalld roplet with high concentration
        V2_0 = 0.00001*V;
        phiA2_0 = 0.45;
        phiB2_0 = 0.45;
        
        S2_0 = (1-phiA2_0-phiB2_0)*V2_0/v;
        A2_0 = phiA2_0*V2_0/(nA*v);
        B2_0 = phiB2_0*V2_0/(nA*v);
        
        [stateOpt, fVal, exitflag] = fmincon(fun, [a-A2_0, b-B2_0, s-S2_0], [], [], [], [], [0, 0, 0], [a, b, s], [], options);
        
        phiA1Opt(k, j) = stateOpt(1) * nA / (stateOpt(1)*nA + stateOpt(2)*nA + stateOpt(3));
        phiB1Opt(k, j) = stateOpt(2) * nA / (stateOpt(1)*nA + stateOpt(2)*nA + stateOpt(3));
        V2Opt(k, j) = V - (stateOpt(1)*nA + stateOpt(2)*nA + stateOpt(3))*v;
        phiA2Opt(k, j) = (a - stateOpt(1))*nA*v / V2Opt(k, j);
        phiB2Opt(k, j) = (b - stateOpt(2))*nA*v / V2Opt(k, j);
        %phiA2Opt(k) = (n - stateOpt(1))*nA*v / V2Opt(k);

        %if (mod(j, 2)==1)
            fprintf('Finished iteration %d %d (phiTotA=%f, phiTotB=%f, exitFlag=%d).\n', k, j, phiTotA, phiTotB, exitflag);
        %end
    end
end

fprintf('done.\n');