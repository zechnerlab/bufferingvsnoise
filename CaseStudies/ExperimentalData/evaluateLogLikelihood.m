function [logL, phi1Pred] = evaluateLogLikelihood(phi, phi1, phiStar, kD_k2_Ratio, varMeanTot, varMeanDilute, covMeanTotMeanDilute)

    %convert parameters in to offset and slope of a linear model
    a = phiStar * kD_k2_Ratio / (1 + kD_k2_Ratio);
    b = 1 / (1 + kD_k2_Ratio);
    
    phi1Pred = gApprox(phi, phiStar, kD_k2_Ratio);
    
    %effective variance needs to be calculated as both x and y axis have
    %uncertainty
    
    effectiveVariance = b^2*varMeanTot + varMeanDilute - 2*b*covMeanTotMeanDilute;
    
    %also need to account for the scaling factor multiplying the
    %exponential in the Gaussian as the standard deviation depends on the
    %parameters too
    logL = sum(-1/2*(phi1 - phi1Pred).^2 ./ (effectiveVariance) - 1/2*log(effectiveVariance));
    
end