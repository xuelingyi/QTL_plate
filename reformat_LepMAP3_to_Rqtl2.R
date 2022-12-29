#### this script was used to reformat the outputs from LepMAP3 for the Rqtl2 analyses ####
### folder p2_215_qtl includes the genotypes (by chr) output from the LepMAP3 script map2genotypes.awk 

## recode genotypes: change 1 to A and 2 to B
#mkdir temp
#for i in {1..21}
#do 
#cat p2_215_qtl/qtldata.LG${i} | sed 's/.*|//g; s/1 1/AA/g; s/1 2/AB/g; s/2 1/BA/g; s/2 2/BB/g' > temp/LG${i}.rqtl
#done

### for better clarity in the plate phenotype mapping: 
# 1 corresponds to the male allele (F1 male is the outbred P. sinensis; replate 1 (A) with S
# 2 corresponds to the female allele (F1 is the hybrid female of P. sinensis and P. pungitius; replate 2 (B) with P)
cat plateQTL_geno_203.csv | sed 's/"AA"/"SS"/g; s/"AB"/"SP"/g; s/"BB"/"PP"/g' > plateQTL_geno_203_SP.csv

#### R script #####
## extract individual names from the pedigree file so that they have the same order as genotypes
pedigree = read.table("./LepMAP3/pedigree_t.txt")
ID = t(pedigree[2,])
ID = ID[-grep("-F0-", ID, fixed = T)]  #no grandparents
ID = ID[-grep("-F1-", ID, fixed = T)]  #no parents
ID = gsub("BC-", "", ID) #match names

### reformat LG1
chr=1

d=read.table(paste0("temp/LG", chr, ".rqtl"), sep="\t")

data = rbind(ID, d[, c(1, 4:ncol(d))])  #remove cols of 0 (chr and male_positions)
data = as.data.frame(t(data))
data[1,1] = "id"
names(data) = data[1,]

nrow(data)
#[1] 217 #number of individuals + marker_row + pos_row

## the reformated genotype data
plateQTL_geno = data[3:217,] # no headers: marker and pos
sites = ncol(data) -1 #colum 1 is ind name

gen = as.data.frame(t(data[1:2,])) #marker and pos
names(gen) = gen[1,] 
gen = gen[-1,]
gen$pos = gen$POS
gen = gen[order(gen$id, decreasing=T),]
names(gen) = c("marker", "chr", "pos")
gen$chr = chr

## the reformated genetic map
## markers were named by their physical positions
plateQTL_gmap = gen


### loop through all the other chr
for (chr in 2:21) {
d=read.table(paste0("temp/LG", chr, ".rqtl"), sep="\t")

data = rbind(ID, d[, c(1, 4:ncol(d))])  #remove cols of 0 (pos and male_pos)
data = as.data.frame(t(data))
data[1,1] = "id"
names(data) = data[1,]

nrow(data)
#[1] 217 #number of individuals
sites = sites + ncol(data) -1 

plateQTL_geno = merge(plateQTL_geno, data[3:217,], by="id", all=T)

gen = as.data.frame(t(data[1:2,])) #marker and pos
names(gen) = gen[1,] 
gen = gen[-1,]
gen$pos = gen$POS
gen = gen[order(gen$id, decreasing=T),]
names(gen) = c("marker", "chr", "pos")
gen$chr = chr
 
plateQTL_gmap = rbind(plateQTL_gmap, gen)

}

## double-check total number of sites
sites

### write out the reformatted data
write.csv(plateQTL_geno, "plateQTL_geno.csv", row.names = F)
write.csv(plateQTL_gmap, "plateQTL_gmap.csv", row.names = F)
