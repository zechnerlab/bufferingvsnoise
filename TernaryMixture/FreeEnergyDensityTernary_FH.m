function f = FreeEnergyDensityTernary_FH(phiA, phiB, chiAA, chiAB, chiBB, kT, v, nA)


    %check if volume fractions are zero, in this case, the phi*log(phi)
    %terms would be nan but should be 0.
    if (phiA == 0)
        sA = 0;
    else
        sA = phiA/nA*log(phiA);    
    end
    
    if (phiB == 0)
        sB = 0;
    else
        sB = phiB/nA*log(phiB);    
    end
    
    if (phiA + phiB == 1)
        sS = 0;
    else
        sS = (1-phiA-phiB)*log(1-phiA-phiB);
    end
    
    
    f = kT/v * (chiAA*phiA*(1-phiA-phiB) + chiAB*phiB*phiA + chiBB*phiB*(1 - phiA-phiB) + sA + sB + sS);

end