function Pi = OsmoticPressureBinary_FH_Approx(phi2, chi, nA)


    Pi = -((phi2*(-1 + nA + nA*phi2*chi) + nA*log(1 - phi2)))/nA;

end