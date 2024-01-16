%%Find minima of binary Flory-Huggins free energy density using fmincon
function [phi1, phi2] = FindBinodalpointsBinary_FH(chi, kT, v, nA, minTol)

    if (nargin<5)
       % This parameter sets the tolerance to call two found minima to be
       % identical. This is important if many multiple different minima
       % close to each other are found.
       minTol = 1e-5; 
    end
    
    options = optimoptions('fmincon', 'Algorithm', 'interior-point');
    options.Display = 'None';
    options.StepTolerance = 1e-12;

    fun = @(phi) FreeEnergyDensityBinary_FH(phi, chi, kT, v, nA);
    
    problem = createOptimProblem('fmincon',...
    'objective',fun,...
    'x0',1e-5,'lb', 0, 'ub', 1, 'options',...
    options);

    ms = MultiStart;
    ms.Display = 'off';
    [x,fval,eflag,output,manymins] = run(ms, problem, 50);

    numMins = length(manymins);
    
    %This calculates the found minima. If multiple minima were found, it uses
    %the minimum and maximum. 
    if (numMins == 1)
        phi1 = manymins(1).X;
        phi2 = phi1;
    elseif (numMins == 2)
        mins = [manymins(1).X, manymins(2).X];
        phi1 = min(mins);
        phi2 = max(mins);
    else
        mins = zeros(size(manymins));
        
        for u=1:length(manymins)
           mins(u) = manymins(u).X; 
        end
        
       mins = uniquetol(mins, minTol);
       
       if (length(mins)>=2)
          phi1 = min(mins);
          phi2 = max(mins);
       elseif (length(mins)==1)
          phi1 = mins;
          phi2 = mins;
       else    
        error('Error happened, more than 2 local minima!\n;');  
       end
    end