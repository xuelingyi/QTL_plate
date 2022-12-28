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
## parameter settings were the same when running on the 215 (complete) and 203 (confirmed sex) indiviudals, except for the maxit variable specified below

### set the inputs
data.name = "sex"  #name of the otuput png and RData files
my.model="binary" 
covar.matrix = NULL
nperm=1000

## the default maxit=100 was used when running on all 215 individuals
## when running on the 203 individuals, set larger maxit to help binary regression to converge
my.maxit=120

print(paste0("nperm: ", nperm, " tol: ", my.tol, " maxit: ", my.maxit))

library(qtl2, lib="../")
plate <- read_cross2("plateQTL.zip")
map <- insert_pseudomarkers(plate$gmap, step=1)
pr <- calc_genoprob(plate, map, error_prob=0.002, cores=0)

m="HK"
  print(paste0("start to run method: ", m))
  # kinship
  k = NULL #default
  if (m=="LMM") {k=kinship}
  if (m=="LOCO") {k=kinship_loco}
  # Performing a genome scan, the output of scan1() is a matrix of LOD scores, positions × phenotypes.
  out <- scan1(genoprobs=pr, pheno=plate$pheno, model=my.model,
               kinship=k, addcovar=covar.matrix, cores=0, tol=my.tol, maxit=my.maxit)
  ymx <- maxlod(out) # overall maximum LOD score
  ymx
  assign(paste0("out_", m), out)
  assign(paste0("ymx_", m), ymx)
  # performing a permutation test
  operm <- scan1perm(genoprobs=pr, pheno=plate$pheno, model=my.model,
                     kinship=k, addcovar=covar.matrix, n_perm=nperm, cores=0, tol=my.tol, maxit=my.maxit)
  summary(operm, alpha=c(0.01, 0.05))
  p01 = summary(operm, alpha=0.01)
  p05 = summary(operm, alpha=0.05)
  assign(paste0("operm_", m), operm)
  assign(paste0("p01_", m), p01)
  assign(paste0("p05_", m), p05)
  # Finding LOD peaks, identify threshold according to previous plots
  # peakdrop indicates the amount that the LOD curve must drop below the lowest of two adjacent peaks
  find_peaks(out, map, threshold=p05[1], peakdrop=1.8, drop=0.1)

  ## plots
pdf(paste0(data.name, ".pdf"))
  y.max = max(ymx, p01[1])
  plot(out, map, col="slateblue", ylim=c(0, y.max*1.35))
  abline(h=p01[1], col="red3")
  abline(h=p05[1], col="black")
  legend("topleft",
         legend=c(colnames(out), "P = 0.01", "P = 0.05"),
         col=c("slateblue", "red3", "black"),
         lwd=c(2,1,1),
         bg="gray90", adj = c(0,0), pt.cex = 0,
         x.intersp = y.max*0.1, y.intersp = y.max*0.1)
dev.off()

## save all data
save(list=c("my.model", "map", "nperm",
            ls(pattern = "out_"), ls(pattern = "ymx_"), ls(pattern = "operm_"), ls(pattern = "p01_"), ls(pattern = "p05_")),
     file=paste0(data.name, ".RData"))
     
     
