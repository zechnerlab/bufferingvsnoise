function f = FreeEnergyDensityBinary_FH(phi, chi, kT, v, nA)


    %check if volume fractions are zero, in this case, the phi*log(phi)
    %terms would be nan but should be 0.
    if (phi == 0)
        sA = 0;
    else
        sA = phi/nA*log(phi);    
    end
    
    if (phi == 1)
        sB = 0;
    else
        sB = (1-phi)*log(1-phi);
    end
    
    
    f = kT/v * (chi*phi*(1-phi) + sA + sB);

end