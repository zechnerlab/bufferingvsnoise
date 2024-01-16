function [phi1Opt, varPhi1, phi2Opt, varPhi2, V2Opt] = CalculateConcentrationDependence_BinaryFH(phiTotGrid, chi, kT, v, nA, V, gamma, output)

if (nargin<8)
   output = 0; 
end


% Initialize vectors for the calculated solutions
phi1Opt = zeros(size(phiTotGrid));
phi2Opt = zeros(size(phiTotGrid));
V2Opt = zeros(size(phiTotGrid));
varPhi1 = zeros(size(phiTotGrid));
varPhi2 = zeros(size(phiTotGrid));

%Intitialize optimizer
options = optimoptions('fmincon', 'Algorithm', 'interior-point');
options.Display = 'None';

if (output==1)
    fprintf('\nCalculating equilibrium concentration-dependence of binary FH model...\n'); 
end

for k=1:length(phiTotGrid)
    % Set total volume fraction in the system and calculate how many total
    % solutes and solvents this corresponds to
    phiTot = phiTotGrid(k);
    n = V*phiTot/(nA*v);
    s = (V - n*nA*v)/v;
    
    % Set function handle to FH Free energy
    fun = @(state) FreeEnergyBinary_FH(state(1), state(2), chi, kT, v, nA, V, n, gamma);

    % initialize a small droplet with high concentration
    V2_0 = 0.0001*V;
    phi2_0 = 0.9;
    
    n1_0 = n - (phi2_0*V2_0)/(nA*v);
    s1_0 = (-n*nA*v + V + (-1 + phi2_0)*V2_0)/v;
    
    % Find minimum of the Free Energy
    stateOpt = fmincon(fun, [n1_0, s1_0], [], [], [], [], [0, 0], [n, s], [], options);
    
%     
%     options = optimoptions('fmincon', 'Algorithm', 'interior-point');
%     options.Display = 'None';
%     options.StepTolerance = 1e-5;
% 
% 
%     problem = createOptimProblem('fmincon',...
%     'objective',fun,...
%     'x0',[n/2, s/2],'lb', [0, 0], 'ub', [n, s], 'options',...
%     options);
% 
%     ms = MultiStart;
%     ms.Display = 'off';
%     [x,fval,eflag,~,manymins] = run(ms, problem, 15);
%     
%     
%     for i=1:length(manymins)
%         state = manymins(i).X;
%         
%         phi1Vec(i) = state(1) * nA / (state(1)*nA + state(2));
%     end
    
    phi1Opt(k) = stateOpt(1) * nA / (stateOpt(1)*nA + stateOpt(2));
    V2Opt(k) = V - (stateOpt(1)*nA + stateOpt(2))*v;
    phi2Opt(k) = (n - stateOpt(1))*nA*v / V2Opt(k);
    
    [~, C_phi] = FreeEnergyBinary_FH(stateOpt(1), stateOpt(2), chi, kT, v, nA, V, n, gamma);
    Cov_Phi = inv(C_phi);
    varPhi1(k) = Cov_Phi(1, 1);
    varPhi2(k) = Cov_Phi(2, 2);
    
    if (mod(k, 10)==1 && output==1)
        fprintf('Finished iteration %d of %d (phiTot=%f).\n', k, length(phiTotGrid), phiTot);
    end
end

if (output==1)
    fprintf('done.\n');
end