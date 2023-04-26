# QTL_sticklebacks_backcross
The R scripts used for the QTL analyses on the backcross Pungitius sticklebacks. "Heterogeneous genomic architecture of skeletal armour traits in sticklebacks".

The linkage maps were constructed using LepMAP3 (Rastas 2017). All QTL mapping was conducted using the R package R/qtl2 (Broman et al 2019). The Lasso regression method was modified from Kemppainen et al. (2021) thanks to the help of Dr. Kemppainen (https://github.com/petrikemppainen).

The zipped Data folder contained input files for these analyses, including: the pedigree, the transformed phenotypes (qtl_pheno_203_final.csv), the formatted genotypes (plateQTL_geno_203_SP.csv), the genetic map (plateQTL_gmap.csv), the physical map (plateQTL_pmap.csv), and the R/qtl2 control file (PlateQTL_203_final.yaml). 


Citations

Broman, K. W., Gatti, D. M., Simecek, P., Furlotte, N. A., Prins, P., Sen, Ś., ... & Churchill, G. A. (2019). R/qtl2: software for mapping quantitative trait loci with high-dimensional data and multiparent populations. Genetics, 211(2), 495-502.

Kemppainen, P., Li, Z., Rastas, P., Löytynoja, A., Fang, B., Yang, J., ... & Merilä, J. (2021). Genetic population structure constrains local adaptation in sticklebacks. Molecular Ecology, 30(9), 1946-1961. 

Rastas, P. (2017). Lep-MAP3: robust linkage mapping even for low-coverage whole genome sequencing data. Bioinformatics, 33(23), 3726-3732.

