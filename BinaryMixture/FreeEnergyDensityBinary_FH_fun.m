function f = FreeEnergyDensityBinary_FH_fun(phi, chi, kT, v, nA)


    %check if volume fractions are zero, in this case, the phi*log(phi)
    %terms would be nan but should be 0.

    sA = phi/nA.*log(phi);    
    sA(phi==0) = 0;
    
    sB = (1-phi).*log(1-phi);
    sB(phi==1) = 0;
    
    
    f = kT/v * (chi*phi.*(1-phi) + sA + sB);

end