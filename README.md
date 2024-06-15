# QTL_sticklebacks_backcross
The R scripts used for the QTL analyses on the backcross Pungitius sticklebacks. Please cite these codes using the Zenodo DOI: 10.5281/zenodo.11559306. The related study "Heterogeneous genomic architecture of skeletal armour traits in sticklebacks" will be published in the Journal of Evolutionary Biology (link to add).

phenotype.R: This script includes the analyses (e.g., correlation) and transformation of phenotypes.

reformat_LepMAP3_to_Rqtl2.R: This script reformat the outputs from LepMAP3 (Rastas 2017) for the Rqtl2 (Broman et al 2019) analyses.

Rqtl2_sex.R: This script is the QTL analyses on the binary trait sex. 

Rqtl2_QTL_mapping.R: This script is the QTL analyses on all the ohter phenotypic traits. 

PVE_LASSO.R: This is the Lasso regression used to esitmate heritabilities. The method was modified from Kemppainen et al. (2021; https://doi.org/10.1111/MEC.15808) and was co-written by Dr. Petri Kemppainen (https://github.com/petrikemppainen).

Input data and files for these analyses can be found in the Dryad archive for this work (link to add).
  
Citations

Broman, K. W., Gatti, D. M., Simecek, P., Furlotte, N. A., Prins, P., Sen, Ś., ... & Churchill, G. A. (2019). R/qtl2: software for mapping quantitative trait loci with high-dimensional data and multiparent populations. Genetics, 211(2), 495-502.

Kemppainen, P., Li, Z., Rastas, P., Löytynoja, A., Fang, B., Yang, J., ... & Merilä, J. (2021). Genetic population structure constrains local adaptation in sticklebacks. Molecular Ecology, 30(9), 1946-1961. 

Rastas, P. (2017). Lep-MAP3: robust linkage mapping even for low-coverage whole genome sequencing data. Bioinformatics, 33(23), 3726-3732.

