# Concentration buffering and noise reduction in phase separating systems
This repository contains Matlab code associated with the manuscript Zechner & Jülicher, "Concentration buffering and noise reduction in non-equilibrium phase separating systems". Note that the title of the manuscript changed during the review-process. The code was used with MATLAB version 9.10.0.1602886 (R2021a) (The MathWorks Inc., Natick, Massachusetts).

## Reproducing the figures
Matlab scripts for reproducing the results in the paper are provided in the folder "CaseStudies" and the respective subfolders. As an example, the folder "CaseStudies/BinaryMixture/" contains all scripts necessary to recreate the results relating to the two-component system. The name of the scripts start with the respective figure number (e.g., "Figure2_*").

## Experimental data
The paper uses previously published concentration measurements of intracellular condensates. Measurements of the 2NTDDX4-YFP system were originally published in Figure 3C and Supplementary Figure S22 in Klosin et al., Science, 2020 (https://doi.org/10.1126/science.aav6691). Measurements of endogenously labelled NPM1 (NPM1-NeonGreen) were published in the same study in Figure 3f. The data has been included in the folder "CaseStudies/ExperimentalData/" with permission from Adam Klosin and Anthony A Hyman. 

Measurements of transiently expressed NPM1-mCherry were originally published in Figure 1b in Riback et al., Nature, 2020 (https://doi.org/10.1038/s41586-020-2256-2). These data are publicly available at: https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2256-2/MediaObjects/41586_2020_2256_MOESM4_ESM.xlsx and a copy of the data has been included in the folder "CaseStudies/ExperimentalData/" (downloaded January 15th 2024).
