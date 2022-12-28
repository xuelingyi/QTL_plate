### QTL mapping on sex (binary trait)
############# prelinimary analyses using all 215 individuals ###############
library(qtl2, lib="./")

#Calculating genotype probabilities
plate <- read_cross2("plateQTL_sex.zip")
nperm=1000
data.name = "qtl_sex"

# insert pseudomarkers into the genetic map
map <- insert_pseudomarkers(plate$gmap, step=1)
#calculate the QTL genotype probabilities,could speed up using multiple cores by adding an argument "cores=x"
pr <- calc_genoprob(plate, map, error_prob=0.002)

# Performing a genome scan, the output of scan1() is a matrix of LOD scores, positions × phenotypes.
out_sex <- scan1(pr, plate$pheno, model="binary")
ymx <- maxlod(out_sex) # overall maximum LOD score
ymx

pdf(file=paste0(data.name, ".pdf"))
plot(out_sex, map, col="slateblue", ylim=c(0, ymx*1.35))#change ymx figure accordingly
legend("topleft", lwd=2, col="slateblue", colnames(out_sex), bg="gray90")
dev.off()

# performing a permutation test
operm <- scan1perm(pr, plate$pheno, model="binary",n_perm=nperm)
summary(operm)
summary(operm, alpha=c(0.01, 0.05))
save(operm, file=paste0(data.name, ".RData"))

# Finding LOD peaks, identify threshold according to previous plots 
# peakdrop indicates the amount that the LOD curve must drop below the lowest of two adjacent peaks
find_peaks(out_sex, map, threshold=0.3, peakdrop=1.8, drop=0.1) 
## outputs:
## lodindex lodcolumn chr    pos      lod  ci_lo  ci_hi
#1       Sex  12 11.286 45.58263 10.819 11.754



############# final: using the 203 individuals ###################
library(qtl2, lib="./")
sex = read_cross2("sex.yaml")
map= insert_pseudomarkers(sex$gmap, step=1)
pr = calc_genoprob(sex, map, error_prob=0.002, cores=0)
out = scan1(genoprobs=pr, pheno=sex$pheno[, "Sex"], model="binary", cores=0, maxit=120)
c2eff = scan1coef(pr[,unique(p$chr)], sex$pheno[, t], contrasts=cbind(mu=c(1,1), a=c(-0.5, 0.5)))
## only dominant effects

# performing a permutation test
operm <- scan1perm(genoprobs=pr, pheno=sex$pheno[, "Sex"], model="binary", cores=0, maxit=120, n_perm=1000)

save(sex, map, pr, out, c2eff, operm,
file=paste0(data, ".RData"))
