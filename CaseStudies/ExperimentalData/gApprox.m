function phi1 = gApprox(phi, phiStar, kD_k2_Ratio)

    %convert parameters in to offset and slope of a linear model
    a = phiStar * kD_k2_Ratio / (1 + kD_k2_Ratio);
    b = 1 / (1 + kD_k2_Ratio);
   
    phi1 = phi;
    phi1(phi>phiStar) = a + b*phi(phi>phiStar);

end