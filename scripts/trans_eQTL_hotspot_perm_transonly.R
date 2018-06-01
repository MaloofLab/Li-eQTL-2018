library(tidyverse)
library(qtl)
library(qtlhot)
library(snowfall) 

load("./output/scanone-eqtl_F2_flipped.RData")
scanone_eQTL.F2 %>% dim() # 4887 56182 

### check correlation, later... 
alphas <- seq(0.01, 0.10, by=0.01)
alphas

### LOD threshold level for e-trait
lod.thr <- lod.thrs[5]
lod.thr 

### trans-eQTL only 
load("./input/cis_trans_result_new_flipped.Rdata") 

cis_eQTL %>% dim() # 11384    12
trans_eQTL %>% dim() # 15197    12  

# include genes with trans-eQTL only  
scanone_eQTL.F2.trans <- scanone_eQTL.F2[,(colnames(scanone_eQTL.F2) %in% c("chr", "pos", trans_eQTL$gene_ID))] 

### get only significant intervals for each e-trait, using LOD drop method 
high1 <- highlod(scanone_eQTL.F2.trans, lod.thr = min(lod.thrs), drop.lod = 1.5)
max(high1, lod.thr = lod.thrs) # max number of e-trait fall into loci with different lod threshold 

hots1 <- hotsize(high1, lod.thr = lod.thr) 
summary(hots1) # for each genomic position, 

# plot(hots1, cex.lab = 1.5, cex.axis = 1.5) 

### permutation to get the statistical significance # permutation takes time, tried to 
set.seed(12345) 

cross.F2$pheno <- cross.F2$pheno[,(colnames(cross.F2$pheno) %in% c(trans_eQTL$gene_ID))]

system.time(
hotperm1.trans <- hotperm(cross = cross.F2, 
n.quant = 600,
n.perm = 100,
lod.thrs = lod.thrs,
alpha.levels = alphas,
drop.lod = 1.5,
verbose = FALSE)
)   

save(hotperm1.trans, file = "./output/hotperm1.trans.Rdata")
