### this script was used for the analyses and transformation of phenotypes 

## original phenotype measures
phenotype = read.csv("plateQTL_phenotype_archive.csv")

## transformed phenotypes for qtl mapping (output from the transform and prep step)
# qtl_pheno = read.csv("qtl_pheno_203_asy.csv") 

library(ggplot2)
library(ggpubr)

################ transform and prep phenotype ########
### take the left-right mean measurements 
phenotype$plateN_mean = (phenotype$PlateN_Left + phenotype$PlateN_Right)/2
phenotype$Area_mean = (phenotype$Left_Plate1_Area + phenotype$Right_Plate1_Area) /2
phenotype$Height_mean = (phenotype$Left_Plate1_Height + phenotype$Right_Plate1_Height) /2
phenotype$Width_mean = (phenotype$Left_Plate1_Width + phenotype$Right_Plate1_Width) /2
phenotype$platedmyomere = (phenotype$PlateMyom.percent_Left + phenotype$PlateMyom.percent_Right) /2

## transform plateN_binary
phenotype$plateN_binary = NULL
for(i in 1:nrow(phenotype)) {
  N = phenotype[i, "plateN_mean"]
  if(N <= 25) {phenotype[i, "plateN_binary"] = "0"} #low-plate
  if(N > 25) {phenotype[i, "plateN_binary"] = "1"} #high-plate
}

phenotype$Plateness = NULL
for (i in 1:nrow(phenotype)) {
  phenotype[i, "Plateness"] = ifelse(phenotype[i, "PlateMyom.percent_Right"] ==100 & 
                                       phenotype[i, "PlateMyom.percent_Left"] ==100, 
                                     "Full", "Partial")
}

## re-code bianry traits
for (i in 1:nrow(phenotype)) {
  phenotype[i, "Sex"] = gsub("Male", "1", phenotype[i, "Sex"])
  phenotype[i, "Sex"] = gsub("Female", "0", phenotype[i, "Sex"])
  phenotype[i, "Sex"] = gsub("unknown", "NA", phenotype[i, "Sex"])
  
  phenotype[i, "Plateness"] = gsub("Full", "1", phenotype[i, "Plateness"])
  phenotype[i, "Plateness"] = gsub("Partial", "0", phenotype[i, "Plateness"])
  
  phenotype[i, "Clutch"] = gsub("C", "", phenotype[i, "Clutch"])
}
unique(phenotype$Sex)
unique(phenotype$Plateness)
phenotype$Clutch = as.numeric(phenotype$Clutch)
phenotype$Sex = as.numeric(phenotype$Sex)
phenotype$Plateness = as.numeric(phenotype$Plateness)
phenotype$plateN_binary = as.numeric(phenotype$plateN_binary)


## sex does not depend on clutch
chisq.test(y=as.factor(phenotype$Sex), x=as.factor(phenotype$Clutch))

## remove the 12 indieividuals of misidentified or unknown sex id
sex_id_remove = read.table("remove_12ind.txt") 

## write out the final phenotype file for the qtl mapping
write.csv(phenotype[!(phenotype$ID %in% sex_id_remove$V1), 
                    c("ID", "Clutch", "Sex",
                        "BodyLength_Stained.fish", "BodyHeight_Stained.fish",
                        "Area_mean", "Height_mean", "Width_mean",
                        "plate_level", "Plateness", "Asymmetry_signed")],
          "qtl_pheno_203.csv", row.names = F)


####################################################
## transformed phenotypes for qtl mapping (output from the transform and prep step)
qtl_pheno = read.csv("qtl_pheno_203_asy.csv") 


###################### dependency on sex or clutch (categorical) ################
## categorical dependent variables 
chisq.test(y=as.factor(qtl_pheno$Sex), x=as.factor(qtl_pheno$Clutch))
chisq.test(y=as.factor(qtl_pheno$Plateness), x=as.factor(qtl_pheno$Sex))
chisq.test(y=as.factor(qtl_pheno$Plateness), x=as.factor(qtl_pheno$Clutch))
chisq.test(y=as.factor(qtl_pheno$plateN_binary), x=as.factor(qtl_pheno$Sex))
chisq.test(y=as.factor(qtl_pheno$plateN_binary), x=as.factor(qtl_pheno$Clutch))


## numeric dependent variables (body size phenotypes)
summary(aov(BodyLength_Stained.fish ~ as.factor(Sex) + as.factor(Clutch), data=qtl_pheno))
summary(aov(BodyHeight_Stained.fish ~ as.factor(Sex) + as.factor(Clutch), data=qtl_pheno))


###################### dependency on sex/clutch/body_size ################
## plate phenotypes dependency on sex or clutch (categorical)
## then control for significant categorical factors (use residuals) to test additional dependency on body size variables

summary(aov(plateN_mean ~ as.factor(Sex) + as.factor(Clutch), data=qtl_pheno))
summary(lm(glm(plateN_mean ~ as.factor(Sex), data = qtl_pheno)$residuals ~ 
             BodyLength_Stained.fish + BodyHeight_Stained.fish,
           data = qtl_pheno))

summary(aov(Area_mean ~ as.factor(Sex) + as.factor(Clutch), data=qtl_pheno))           
summary(lm(glm(Area_mean ~ as.factor(Sex), data = qtl_pheno)$residuals ~ 
             BodyLength_Stained.fish + BodyHeight_Stained.fish,
           data = qtl_pheno))

summary(aov(Height_mean ~ as.factor(Sex) + as.factor(Clutch), data=qtl_pheno))
summary(lm(Height_mean ~ 
           BodyLength_Stained.fish + BodyHeight_Stained.fish,
           data = qtl_pheno))

summary(aov(Width_mean ~ as.factor(Sex) + as.factor(Clutch), data=qtl_pheno))
summary(lm(glm(Width_mean ~ as.factor(Sex) + as.factor(Clutch), data = qtl_pheno)$residuals ~ 
             BodyLength_Stained.fish + BodyHeight_Stained.fish,
           data = qtl_pheno))

summary(aov(Asymmetry_signed ~ as.factor(Sex) + as.factor(Clutch), data=qtl_pheno))
summary(lm(Asymmetry_signed ~ BodyLength_Stained.fish + BodyHeight_Stained.fish,
           data = qtl_pheno))

############# check distribution #############
png("versions/figure/distribution_final.png", width = 8, height = 9, units = "in", res=600)
pdf("distribution_final.pdf", width = 8, height = 9)
par(mfrow=c(4,4))
for (trait in c("Length", "Height")) {
  my.data = phenotype[, paste0("Body", trait, "_Stained.fish")]
  hist(my.data, main = paste0("Body ", trait),
       xlab=NULL)
  qqnorm(my.data, pch = 1, frame = FALSE, xlab=paste0("Body ", trait))
  qqline(my.data, col = "steelblue", lwd = 2)
}
hist(phenotype$plateN_mean, main = "plate number mean",
     xlab=NULL)
qqnorm(phenotype$plateN_mean, pch = 1, frame = FALSE,xlab="plate number mean")
qqline(phenotype$plateN_mean, col = "steelblue", lwd = 2)
hist(phenotype$platedmyomere, main = "mean proportion of\nplated myomere",
     xlab=NULL)
qqnorm(phenotype$platedmyomere, pch = 1, frame = FALSE,xlab="mean proportion of\nplated myomere")
qqline(phenotype$platedmyomere, col = "steelblue", lwd = 2)
hist(phenotype$plate_level, main = "plate level", xlab=NULL)
qqnorm(phenotype$plate_level, pch = 1, frame = FALSE,xlab="plate level")
qqline(phenotype$plate_level, col = "steelblue", lwd = 2)
for (trait in c("Area", "Height", "Width")) {
  hist(phenotype[, paste0(trait, "_mean")], 
       xlab=NULL, main=paste0("Plate1_", trait, "_mean"))
  qqnorm(phenotype[, paste0(trait, "_mean")], pch = 1, frame = FALSE,
         xlab=paste0("Plate1_", trait, "_mean"))
  qqline(phenotype[, paste0(trait, "_mean")], col = "steelblue", lwd = 2)
}
dev.off()
par()

######################### sexual dimorphism #################################
pheno203 = read.csv("qtl_pheno_203_final.csv")
pheno203$Sex = as.factor(pheno203$Sex)
library(plotrix)

## sexual dimorphism in body size
t.test(pheno203[pheno203$Sex == "0", "BodyLength_Stained.fish"],
       pheno203[pheno203$Sex == "1", "BodyLength_Stained.fish"])
t.test(pheno203[pheno203$Sex == "0", "BodyHeight_Stained.fish"],
       pheno203[pheno203$Sex == "1", "BodyHeight_Stained.fish"])

## sexual dimorphism in plate number and size
t.test(pheno203[pheno203$Sex == "1", "plateN_mean"], 
       pheno203[pheno203$Sex == "0", "plateN_mean"])
t.test(pheno203[pheno203$Sex == "1", "Area_mean"], 
       pheno203[pheno203$Sex == "0", "Area_mean"])

## control for body length
pheno203$residual = lm(plateN_mean ~ BodyLength_Stained.fish, data = pheno203)$residuals
t.test(pheno203[pheno203$Sex == "1", "residual"], 
       pheno203[pheno203$Sex == "0", "residual"])

pheno203$residual = lm(Area_mean ~ BodyLength_Stained.fish, data = pheno203)$residuals
t.test(pheno203[pheno203$Sex == "1", "residual"], 
       pheno203[pheno203$Sex == "0", "residual"])








