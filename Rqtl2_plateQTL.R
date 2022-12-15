####### notes
## one folder was created for each QTL analysis; folders were named by the QTL analyses (e.g., height)
## each folder contains one input file plateQTL.zip; each plateQTL.zip file has the same gmap data, the same geno data, but different phenotypes, covariates, and control files corresponding to each QTL analyses
## the following R scripts were run in the same directory as the QTL folders

####### R scripts 
## set the inputs: modify accordingly for each QTL analyses (example below is the height QTL)
## all QTL model settings are available in the corresponding publication 

data.name = "height"  # name of the QTL analyses (same as folder names)
# all data.names=c("plateness", "length", "height", "platepercent",  "plate1area", "asymmetry")

my.model="normal" #for numeric variables, least squares
#my.model="binary" #for binary variables, logistic regression

my.covar = c("Sex", "BodyLength_Stained.fish")
#my.covar = NULL ## if no covariates
# covariates choose from "Sex", "BodyLength_Stained.fish", "BodyHeight_Stained.fish", "Plate1_Area_mean", "Asymmetry", "Plateness", "PlateMyom_percent_mean"

my.method = c("HK", "LMM", "LOCO")
#my.method = c("HK") # if binary models

nperm=1000  #number of permutation


#### no modification needed in the scripts below
# double-check my input settings  
print(paste0("data.name: ", data.name))
print(paste0("my.model: ", my.model))
print(paste0("my.covar: ", paste(my.covar, collapse = ", ")))
print(paste0("nperm: ", nperm))

library(qtl2)
## go to the QTL folder
setwd(data.name)

# load the data
plate <- read_cross2("plateQTL.zip")
# insert pseudomarkers into the genetic map
map <- insert_pseudomarkers(plate$gmap, step=1)
#Calculating genotype probabilities
pr <- calc_genoprob(plate, map, error_prob=0.002)

## kinship matrix for LMM and LOCO methods
#Calculate a kinship matrix
kinship <- calc_kinship(pr)
# use the “leave one chromosome out” (LOCO) method scan each chromosome using a kinship matrix that is calculated using data from all other chromosomes
kinship_loco <- calc_kinship(pr, "loco")

## get covariants
covar.matrix = NULL
if (!(is.null(my.covar))) {
  for (i in my.covar){
    j = which(names(plate$covar) == i) #the column of that variable
    covar.matrix = cbind(covar.matrix,
                         setNames(as.numeric(plate$covar[,j]), rownames(plate$covar)))
  }
}

## loop through my.method
for (m in my.method){
  print(paste0("start to run method: ", m))
  
  # kinship matrix 
  k = NULL #default
  if (m=="LMM") {k=kinship}
  if (m=="LOCO") {k=kinship_loco}
  
  # Performing a genome scan, the output of scan1() is a matrix of LOD scores, positions × phenotypes.
  out <- scan1(genoprobs=pr, pheno=plate$pheno, model=my.model,
               kinship=k, addcovar=covar.matrix)
  ymx <- maxlod(out) # overall maximum LOD score
  print(ymx)
  assign(paste0("out_", m), out)
  assign(paste0("ymx_", m), ymx)
  
  # performing a permutation test
  operm <- scan1perm(genoprobs=pr, pheno=plate$pheno, model=my.model,
                     kinship=k, addcovar=covar.matrix, n_perm=nperm)
  print(summary(operm, alpha=c(0.01, 0.05)))
  p01 = summary(operm, alpha=0.01)
  p05 = summary(operm, alpha=0.05)
  assign(paste0("operm_", m), operm)
  assign(paste0("p01_", m), p01)
  assign(paste0("p05_", m), p05)
  
  # Finding LOD peaks, identify threshold according to previous plots
  # peakdrop indicates the amount that the LOD curve must drop below the lowest of two adjacent peaks
  print(find_peaks(out, map, threshold=p05[1], peakdrop=1.8, drop=0.1))

  ## plots
  y.max = max(ymx, p01[1])
  plot(out, map, col="slateblue", ylim=c(0, y.max*1.35))
  abline(h=p01[1], col="red3")
  abline(h=p05[1], col="black")
  legend("topleft",
         legend=c(colnames(out), "P = 0.01", "P = 0.05"),
         col=c("slateblue", "red3", "black"),
         lwd=c(2,1,1),
         bg="gray90", adj = c(0,0), pt.cex = 0,
         x.intersp = 1, y.intersp = 1)
}
## save all data
save(list=c("my.model", "my.covar", "map", "nperm",
            ls(pattern = "out_"), ls(pattern = "ymx_"), ls(pattern = "operm_"), ls(pattern = "p01_"), ls(pattern = "p05_")),
     file=paste0(data.name, ".RData"))
