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

############# plot results ################
## lateral plates
library(qtl2)
pdf("versions/figure/plate_QTL.pdf", height = 10, width = 8)
png("versions/figure/plate_QTL.png", height = 10, width = 8, units = "in", res = 600)

par(mfrow=c(2,1))

number.traits = c("plateN_mean", "Plateness")
number.colors = c("#113b69", "#2b83ba")
for (i in 1:2) {
  trait = number.traits[i]
  mycolor = ifelse(i>1, alpha(number.colors[i], 0.2), number.colors[i])
  load(paste0("Rqtl/", trait, ".RData"))
  ##[1]data, map, pr, covar.matrix, out, operm
  p05 = summary(operm, alpha=0.05)
  peak = find_peaks(out, map, threshold=p05[1], peakdrop=2, prob = 0.95)
  peak$lodcolumn = trait
  plot(out, map, col=mycolor, ylim=c(0, max(maxlod(out), p05[1])*1.35), 
       add=ifelse(i>1, T, F), type="b", pch=16, cex=0.6, xlab="Linkage Group")
  abline(h=p05[1], col=number.colors[i])}
legend("topleft", 
       legend=c("mean plate number", "Plateness", "p-value = 0.05"),
       col=c(number.colors, "black"),
       pch=c(16,16, 16), pt.cex = c(0.8, 0.8, 0), 
       lwd=c(-0.1, -0.1, 1.5),
       bg="gray95", adj = c(0,0), 
       x.intersp = 1, y.intersp = 1)  
title("A) plate number QTL", adj=0)

size.traits = c("Height_mean", "Area_mean", "Width_mean")
size.colors = c("#192b16", "#6e8f69", "#17ad00")
for (i in 1:3) {
  trait = size.traits[i]
  mycolor = ifelse(i>1, alpha(size.colors[i], 0.2), size.colors[i])
  load(paste0("Rqtl/", trait, ".RData"))
  ##[1]data, map, pr, covar.matrix, out, operm
  p05 = summary(operm, alpha=0.05)
  peak = find_peaks(out, map, threshold=p05[1], peakdrop=2, prob = 0.95)
  peak$lodcolumn = trait
  plot(out, map, col=mycolor, ylim=c(0, max(maxlod(out), p05[1])*1.35), 
       add=ifelse(i>1, T, F), type="b", pch=16, cex=0.6, xlab="Linkage Group")
  abline(h=p05[1], col=size.colors[i])}
legend("topleft", 
       legend=c(gsub("_", " ", size.traits), "p-value = 0.05"),
       col=c(size.colors, "black"),
       pch=c(16,16,16,16), pt.cex = c(0.8, 0.8, 0.8, 0), 
       lwd=c(-0.1, -0.1, -0.1, 1.5),
       bg="gray95", adj = c(0,0), 
       x.intersp = 1, y.intersp = 1)  
title("B) plate size QTL", adj=0)

dev.off()


F2_order = row.names(data$pheno)
pheno = pheno203[pheno203$ID %in% F2_order, ]
pheno = pheno203[order(pheno203$ID, levels=F2_order),]

SNPcol <- c("slateblue", "violetred")
trait = "plateN_mean"
model = "normal"
plateN_mean.eps = pheno[, c("ID", "plateN_mean")]

pdf(paste0("figure/", trait, "_genotypes.pdf"), width = 5, height = 7)
png(paste0("figure/", trait, "_genotypes.png"), width = 5, height = 7, units = "in", res=600)
par(mfrow = c(3,2))
for (chr in unique(peak$chr)) {
  pos = peak[peak$chr == chr, "pos"]
  c2eff = scan1coef(pr[, chr], data$pheno[, trait], model=model)
  plot_coef(c2eff, map[chr], columns = 1:2, col=SNPcol, 
            main=paste0("LG", chr, " pos_", pos, " (cM)"))
  last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients
  for(i in seq(along=last_coef)) {
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=SNPcol[i])}
  
  g = maxmarg(pr, map, chr=chr, pos=pos, return_char = T)
  plot_pxg(g, data$pheno[, trait], ylab=trait, sort=F, pch=16, 
           col=c("violetred", "slateblue")[factor(g)])
  
  names = names(plateN_mean.eps)
  plateN_mean.eps = cbind(plateN_mean.eps, g)
  names(plateN_mean.eps) = c(names, paste0("LG", chr, "_pos", pos))
}
dev.off()
par()


plateN_mean.eps = read.csv("plateN_mean.eps.csv")
pdf("versions/figure/plateN_mean_epistatic_violin2.pdf", width = 10, height = 7)
png("versions/figure/plateN_mean_epistatic_violin2.png", width = 10, height = 7, units="in", res=600)
ggplot(plateN_mean.eps) + 
  geom_violin(aes(x=Genotype, y=plateN_mean), scale="width") +
  geom_jitter(aes(x=Genotype, y=plateN_mean), width=0.05) +
  scale_x_discrete(limits=c("SP/SP/SS", "SP/SS/SS", "SS/SP/SS", "SP/SP/SP", 
                            "SS/SS/SS", "SS/SP/SP", "SP/SS/SP", "SS/SS/SP")) +
  xlab("Genotype on LG4/LG13/LG17 QTL peaks") +
  ylab("Plate number")
dev.off()


### body size
library(qtl2)
my.traits = c("BodyLength_Stained.fish", "BodyHeight_Stained.fish")
body.colors = c("#311800", "#f48a21")

pdf("versions/figure/bodysize_QTL_v2LG.pdf", height = 5, width = 8)
png("versions/figure/bodysize_QTL_v2LG.png", height = 5, width = 8, units = "in", res = 600)

for (i in 1:2) {
  trait = my.traits[i]
  mycolor = ifelse(i>1, alpha(body.colors[i], 0.3), body.colors[i])
  
  load(paste0("Rqtl/", trait, (".RData")))
  p05 = summary(operm, alpha=0.05)
  peak = find_peaks(out, map, threshold=p05[1], peakdrop=2, prob = 0.95)
  peak$lodcolumn = trait
  
  plot(out, map, col=mycolor, ylim=c(0, max(maxlod(out), p05[1])*1.35), 
       add=ifelse(i>1, T, F), type="b", pch=16, cex=0.6,
       xlab="Linkage Group")
       #xlab="", ylab="")
  abline(h=p05[1], col=body.colors[i])
}
legend("topleft", 
       legend=c("Standard Length", "Body Height", "p-value = 0.05"),
       col=c(body.colors, "black"),
       pch=c(16,16, 16), pt.cex = c(0.8, 0.8, 0), 
       lwd=c(-0.1, -0.1, 1.5),
       bg="gray95", adj = c(0,0), 
       x.intersp = 1, y.intersp = 1)  
dev.off()

## all maps
traits = c("BodyLength_Stained.fish", "BodyHeight_Stained.fish", 
           "plateN_mean", "plateN_binary25", 
           "Plateness", 
           "Area_mean", "Height_mean", "Width_mean")

pdf("versions/figure/QTL_maps_final.pdf", height = 14, width = 12)
png("versions/figure/QTL_maps_final.png", height = 14, width = 12, units = "in", res = 600)
par(mfrow=c(4,2))
for (my.trait in traits) {
  load(paste0("Rqtl/", my.trait, (".RData")))
  ##[1]data, map, pr, covar.matrix, out, operm
  p05 = summary(operm, alpha=0.05)
  peak = find_peaks(out, map, threshold=p05[1], peakdrop=2, prob = 0.95)
  plot(out, map, col="slateblue", ylim=c(0, max(maxlod(out), p05[1])*1.35)) 
  abline(h=p05[1], col="black")
  legend("topleft", 
         legend=c(my.trait, "p-value = 0.05"),
         col=c("slateblue", "black"),
         lwd=c(2,1), 
         bg="gray90", adj = c(0,0), pt.cex = 0, 
         x.intersp = 1, y.intersp = 1)  
}
dev.off()
