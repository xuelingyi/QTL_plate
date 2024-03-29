pheno203 = read.csv("qtl_pheno_203_final.csv")

## prep snps
pmap = read.csv("plateQTL_pmap.csv") # 59629 SNPs
pmap$pos = as.numeric(format(pmap$pos * 1000000, scientific=F)) # change the unit from megabase to bp

geno203 = read.csv("plateQTL_geno_203_SP.csv")
genoorder = as.data.frame(names(geno203))
# genoorder = read.table("geno203_order.list", header=T)
genoorder = genoorder[-(1),]

## make sure snps are in the same order in pamp (rows) and geno (cols)
pmap$marker = factor(pmap$marker, levels = genoorder)
pmap = pmap[order(pmap$marker), ]
all(pmap$marker == genoorder)
#TRUE

mylist = list.files("/QTL_range/")
#[1] "Area_mean_19_36.939.list"               "Area_mean_20_53.44.list"               
#[3] "BodyHeight_Stained.fish_6_55.305.list"  
#[4] "BodyLength_Stained.fish_12_87.896.list"  "BodyLength_Stained.fish_16_43.088.list" 
#[6] "Height_mean_20_53.44.list"             
#[7] "plateN_mean_QTLrange_LG13.list"  "plateN_mean_QTLrange_LG17.list"    "plateN_mean_QTLrange_LG4.list"          
#[10] "Plateness_13_63.24.list"               
#[11] "sex203_12_11.286.list"                  "sex203_12_13.623.list"    
#[13] "Width_mean_17_65.635.list" "Width_mean_19_33.707.list" "Width_mean_20_53.44.list"      

mytraits = names(pheno203)[c(4:8,10,12)] #no binary traits (sex, plateN_binary, Plateness)

for (trait in mytraits) {
  # code nonQTL regions as 0
  pmap$QTL=0
  # grep all QTL regions (95% CI) of this trait (all traits have QTL on different chr except sex)
  #files = ifelse(trait == "Sex", mylist[grep("sex203", mylist)], mylist[grep(trait, mylist)])
  files = mylist[grep(trait, mylist)]
    
  for (i in files) {
    data = read.table(paste0("./PVELG/QTL_range/", i))
    names(data) = c("LG", "pos_genetic")
    data$LG = gsub("LG", "", data$LG)
    
    if(all(data$pos_genetic %in% pmap[pmap$chr == unique(data$LG), "pos"])) {
      pmap[pmap$chr == unique(data$LG) & pmap$pos %in% data$pos_genetic, "QTL"] = unique(data$LG)
    }else {
      print(paste0(i, ": not all pos in map"))
      break}
    #QTLn = QTLn +1
  }
  print(paste0(trait, ": ", paste0(unique(pmap$QTL), collapse = ", ")))
  assign(paste0(trait, "_QTL"), pmap$QTL)
}
#[1] "BodyLength_Stained.fish: 0, 12, 16"
#[1] "BodyHeight_Stained.fish: 0, 6"
#[1] "Area_mean: 0, 19, 20"
#[1] "Height_mean: 0, 20"
#[1] "Width_mean: 0, 17, 19, 20"
#[1] "Plateness: 0, 13"
#[1] "Asymmetry_signed: 0"
#[1] "plateN_mean: 0, 4, 13, 17"
pmap_QTL = pmap
pmap_QTL = cbind(pmap_QTL, Asymmetry_signed_QTL, Area_mean_QTL, Height_mean_QTL, Width_mean_QTL,
                 BodyHeight_Stained.fish_QTL, BodyLength_Stained.fish_QTL, 
                 plateN_mean_QTL, Plateness_QTL)
write.csv(pmap_QTL, "pmap_QTL.csv", row.names = F)


library(glmnet)
#########################
## codes of the PVELG function was originally provided in Kemppainen, P., Li, Z., Rastas, P., Löytynoja, A., Fang, B., Yang, J., ... & Merilä, J. (2021). Genetic population structure constrains local adaptation in sticklebacks. Molecular Ecology, 30(9), 1946-1961. 
## the PVELG function was slighly modified here to adapt to our dataset
#########################
PVELG <- function(X,y,chromosome){
  #evaluating the narrow sense heritability of the quantitative trait and the percentage of variance (PVE) explained by each  linkage group using the LASSO
  
  #X <- as.matrix(X)
  X <- data.matrix(X)
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  LAM <- NULL
  i <- 1
  lam0 <-mean(sapply(1:15, function(i){
    result <- cv.glmnet(X,y,alpha=1)
    LAM[i] <- result$lambda.min
  }))
  
  gfit <- glmnet(X,y,lambda=lam0,alpha=1)
  htheta  <- gfit$beta
  htheta <- as.numeric(htheta)
  vg <- sum((y-gfit$a0-X%*%htheta)^2)/(n-gfit$df-1)
  
  PVE <-(var(y)-vg)/var(y)
  
  PVE.chro <- NULL
  #PVE.chro2 <- NULL
  
  ID <- as.numeric(unique(chromosome))
  
  for (i in 1:length(ID)){
    vg1 <- sum((y-gfit$a0-as.matrix(X[,-which(chromosome==ID[i])])%*%htheta[-which(chromosome==ID[i])])^2)/(n-sum(abs(htheta[-which(chromosome==ID[i])])>0)-1)
    PVE.chro[i] <- PVE-(var(y)-vg1)/var(y)
    #PVE.chro2[i] <- (var(y)-vg1)/var(y)
  }
  return(list(PVE=PVE, LG=ID, PVE.chro=PVE.chro))
  
} 

## input files
pmap_QTL = read.csv("pmap_QTL.csv")
pheno203 = read.csv("qtl_pheno_203_final.csv")
geno203 = read.csv("plateQTL_geno_203_SP.csv")

## recode genotypes for LASSO
geno203 = as.data.frame(apply(geno203, 2, function(x){gsub("SS", "-1", x)}))
geno203 = as.data.frame(apply(geno203, 2, function(x){gsub("SP", "1", x)}))
geno203 = apply(geno203[, -1], 2, function(x){as.numeric(x)})
dim(geno203)
#[1]   203 59629
# number of individuals, number of all SNPs

## remove duplicated SNPs (loci having the same genotypes)
notDuplicated <- !duplicated(sapply(1:ncol(geno203), function(x){
  paste(geno203[,x],collapse = "")
}))
geno203 = geno203[,notDuplicated]
dim(geno203)
#[1]  203 2677 (nonduplicated SNPs)
dim(pmap_QTL)
#[1] 59629    12
pmap_QTL = pmap_QTL[notDuplicated,]
dim(pmap_QTL)
# [1] 2677   12
all(colnames(geno203) == pmap_QTL$marker)
#[1] TRUE
save(geno203, pmap_QTL, pheno203, file="PVELG_nodup.RData")


load("PVELG_nodup.RData")
mytraits = c("BodyLength_Stained.fish", "BodyHeight_Stained.fish",
             "Area_mean", "Height_mean", "Width_mean",
             "Plateness", "plateN_mean", "Asymmetry_signed")
myPVELG = as.data.frame(NULL)

for (trait in mytraits){
  geno=as.matrix(geno203)
  chr=as.vector(unlist(pmap_QTL[, paste0(trait, "_QTL")]))
  ## sex was not controlled only in Asymmetry and Height_mean
  if(trait %in% c("Asymmetry_signed", "Height_mean")){
    phe = pheno203[, trait]
  } else {phe = resid(lm(pheno203[, trait] ~ pheno203$Sex))}
  phe=as.vector(unlist(phe))
  
  PVE = PVELG(X=geno, y=phe, chromosome=chr)
  print(PVE)
  assign(paste0(trait, "_PVE"), PVE)
  
  output = cbind(as.data.frame(PVE$LG),
                 as.data.frame(round(PVE$PVE.chro,4)))
  names(output) = c("LG", "PVE.chro")
  output$trait = paste0(trait, "::PVE total/sum=", 
                        round(PVE$PVE,4), "/",  
                        round(sum(PVE$PVE.chro),4))
  myPVELG = rbind(myPVELG, output)
}
write.csv(myPVELG, "PVELG_nodup_output.csv", row.names = F)
save(list=ls(pattern = "_PVE"), file="PVELG_nodup_outputs.RData")
