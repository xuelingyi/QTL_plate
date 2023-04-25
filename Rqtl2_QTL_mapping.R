#### These are the r scripts used to conduct QTL mapping on each phenotypic trait (except for sex)

####################### binary traits ####################
my.trait = "Plateness"
# or my.trait = "plateN_binary"

library(qtl2, lib="./")
data=read_cross2("PlateQTL_203_final.yaml")
map=insert_pseudomarkers(data$gmap, step=1)
pr=calc_genoprob(data, map, error_prob=0.002, cores=0)
out=scan1(genoprobs=pr, pheno=data$pheno[, my.trait], model="binary", cores=0)
#performing a permutation test
operm = scan1perm(genoprobs=pr, pheno=data$pheno[, my.trait], model="binary", cores=0, n_perm=1000)
save(data, map, pr, out, operm, file=paste0(my.trait, ".RData"))


####################### the other numeric traits ##############################
### input settings (each trait was run independently):
my.trait = "Area_mean"
my.covar = c("Sex", "BodyLength_Stained.fish")

my.trait = "BodyHeight_Stained.fish"
my.covar = c("Sex", "Clutch")

my.trait = "BodyLength_Stained.fish"
my.covar = c("Sex", "Clutch")

my.trait = "Height_mean"
my.covar = c("BodyLength_Stained.fish")

my.trait = "plateN_mean"
my.covar = c("Sex")

my.trait = "Plateness"

my.trait = "Width_mean"
my.covar = c("Sex", "Clutch", "BodyLength_Stained.fish")

###### source the same codes below
library(qtl2)
data=read_cross2("PlateQTL_203_asy.yaml")
map=insert_pseudomarkers(data$gmap, step=1)
pr=calc_genoprob(data, map, error_prob=0.002, cores=0)

covar.matrix = NULL
if (!(is.null(my.covar))) {
  for (i in my.covar){
    j = which(names(data$covar) == i) #the column of that variable
    covar.matrix = cbind(covar.matrix,
                         setNames(as.numeric(data$covar[,j]), rownames(data$covar)))
  }
}

out=scan1(genoprobs=pr, pheno=data$pheno[, my.trait], model="normal", cores=0, addcovar=covar.matrix)
#performing a permutation test
operm = scan1perm(genoprobs=pr, pheno=data$pheno[, my.trait], model="normal", cores=0, n_perm=1000, addcovar=covar.matrix)
save(data, map, pr, covar.matrix, out, operm, file=paste0(my.trait, ".RData"))
