---
title: "QTL_eQTL_mapping_2"
author: "Ruijuan Li"
date: "5/13/2018"
output: 
  html_document: 
    keep_md: yes
---


```r
library(tidyverse)
```

```
## ── Attaching packages ────────────────────────────────── tidyverse 1.2.1 ──
```

```
## ✔ ggplot2 2.2.1     ✔ purrr   0.2.4
## ✔ tibble  1.4.2     ✔ dplyr   0.7.4
## ✔ tidyr   0.8.0     ✔ stringr 1.3.0
## ✔ readr   1.1.1     ✔ forcats 0.3.0
```

```
## Warning: package 'ggplot2' was built under R version 3.2.5
```

```
## Warning: package 'readr' was built under R version 3.2.5
```

```
## Warning: package 'purrr' was built under R version 3.2.5
```

```
## Warning: package 'dplyr' was built under R version 3.2.5
```

```
## ── Conflicts ───────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
library(qtl) 
```

```
## Warning: package 'qtl' was built under R version 3.2.5
```

```r
library(Biostrings)  
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     grep, grepl, intersect, is.unsorted, lapply, lengths, Map,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unlist, unsplit
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:dplyr':
## 
##     rename
```

```
## Loading required package: IRanges
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
```

```
## The following objects are masked from 'package:purrr':
## 
##     reduce, simplify
```

```
## The following object is masked from 'package:tidyr':
## 
##     expand
```

```
## Loading required package: XVector
```

```
## 
## Attaching package: 'XVector'
```

```
## The following object is masked from 'package:purrr':
## 
##     compact
```

### QTL mapping result output 

```r
setwd("~/Desktop/F2_paper/submission/Li-eQTL-TAG-2018/")
load("input/QTL_result_all.Rdata")

### summarize QTL mapping result  
threshold.95 <- tibble(perm.threshold = bind_rows(scanone.perm.imp.all) %>% as.numeric(), 
                       trait = colnames(bind_rows(scanone.perm.imp.all)))

scanone.qtl.2 <-  
bind_cols(scanone.imp.all) %>% 
  dplyr::select(chr, pos, starts_with("lod"))
rownames(scanone.qtl.2) <- rownames(scanone.imp.all$Crude_oil_contents)
colnames(scanone.qtl.2)[3:ncol(scanone.qtl.2)] <- names(scanone.imp.all)

scanone.gather <- scanone.qtl.2 %>%
  gather(key = trait, value = LOD, -chr, -pos) %>%
  left_join(threshold.95)
```

```
## Joining, by = "trait"
```

```r
sig.chrs <- scanone.gather %>% dplyr::filter(LOD > perm.threshold) %>%
  group_by(trait,chr) %>%
  dplyr::summarise(count = n())

# now for each significant chromosome/trait combo run bayesint
bayesint.list <- apply(sig.chrs,1,function(hit) { # for every row("trait, chr, count") in eigengene module 
    result <- bayesint(scanone.qtl.2[c("chr","pos",hit["trait"])],  
                     chr=hit["chr"], 
                     lodcolumn = 1, 
                     expandtomarkers = TRUE 
  )
  colnames(result)[3] <- "LOD" 
  result
})  
 
names(bayesint.list) <- sig.chrs$trait

bayesint.list <- lapply(bayesint.list,function(x)  
                          x %>% 
                          as.data.frame() %>%
                          rownames_to_column(var="markername")  %>% # make rownames to column and use "markername" as the colname for the new colomn  
                          mutate(chr=as.character(chr))
) 

bayesint.list %>% length() # 33
```

```
## [1] 33
```

```r
bayesint.list.scanone <- bayesint.list 

bayesint.result.scanone <- as.tibble(bind_rows(bayesint.list,.id="trait")) %>% # combine list into tibble 
  dplyr::select(trait,chr,pos,markername,LOD) %>% 
  separate(markername,into=c("chr1","Mbp"),sep="_", convert=TRUE) %>% 
  group_by(trait,chr) %>% 
  dplyr::summarize(start=min(Mbp, na.rm = T),end=max(Mbp, na.rm = T),min_eQTL_LOD=min(LOD),max_eQTL_LOD=max(LOD), genetic_start=min(pos, na.rm = T), genetic_end=max(pos, na.rm = T)) %>% 
  #for the high QTL peaks the interval width is 0.  That is overly precise and need to widen those.
  mutate(start=ifelse(start==end,max(0,start-20000),start), end=ifelse(start==end,end+20000,end))
```

```
## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 15 rows [5,
## 11, 26, 32, 35, 47, 53, 56, 62, 68, 83, 86, 92, 95, 98].
```

```r
bayesint.result.scanone %>% dim() # 33 8 
```

```
## [1] 33  8
```

```r
### cim
threshold.95 <- tibble(perm.threshold = bind_rows(cim.perm.all) %>% as.numeric(), 
                       trait = colnames(bind_rows(cim.perm.all)))

cim.qtl.2 <-  
bind_cols(cim.qtl.all) %>% 
  dplyr::select(chr, pos, starts_with("lod"))

rownames(cim.qtl.2) <- rownames(cim.qtl.all$Crude_oil_contents)
colnames(cim.qtl.2)[3:ncol(cim.qtl.2)] <- names(cim.qtl.all)

cim.gather <- cim.qtl.2 %>%
  gather(key = trait, value = LOD, -chr, -pos) %>%
  left_join(threshold.95)
```

```
## Joining, by = "trait"
```

```r
# look for overlap, for each trait, find QTL border and look for genes under QTL peaks 
sig.chrs <- cim.gather %>% dplyr::filter(LOD > perm.threshold) %>%
  group_by(trait,chr) %>%
  dplyr::summarise(count = n()) 

# now for each significant chromosome/trait combo run bayesint
bayesint.list <- apply(sig.chrs,1,function(hit) { # for every row("trait, chr, count") in eigengene module 
    result <- bayesint(cim.qtl.2[c("chr","pos",hit["trait"])],  
                     chr=hit["chr"], 
                     lodcolumn = 1, 
                     expandtomarkers = TRUE 
  )
  colnames(result)[3] <- "LOD" 
  result
})  

names(bayesint.list) <- sig.chrs$trait

bayesint.list <- lapply(bayesint.list,function(x)  
                          x %>% 
                          as.data.frame() %>%
                          rownames_to_column(var="markername")  %>% # make rownames to column and use "markername" as the colname for the new colomn  
                          mutate(chr=as.character(chr))
) 

bayesint.list %>% length() # 26 
```

```
## [1] 26
```

```r
# save bayesint result for later 
bayesint.list.cim <- bayesint.list

bayesint.result.cim <- as.tibble(bind_rows(bayesint.list,.id="trait")) %>% # combine list into tibble 
  dplyr::select(trait,chr,pos,markername,LOD) %>% 
  separate(markername,into=c("chr1","Mbp"),sep="_", convert=TRUE) %>% 
  group_by(trait,chr) %>% 
  dplyr::summarize(start=min(Mbp, na.rm = T),end=max(Mbp, na.rm = T),min_eQTL_LOD=min(LOD),max_eQTL_LOD=max(LOD), genetic_start=min(pos, na.rm = T), genetic_end=max(pos, na.rm = T)) %>% 
  #for the high QTL peaks the interval width is 0.  That is overly precise and need to widen those.
  mutate(start=ifelse(start==end,max(0,start-20000),start), end=ifelse(start==end,end+20000,end))
```

```
## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 8 rows [2,
## 5, 29, 32, 44, 59, 68, 71].
```

```r
bayesint.result.cim %>% dim() # 26 8
```

```
## [1] 26  8
```

```r
bayesint.result.cim$model <- rep("cim", nrow(bayesint.result.cim))
bayesint.result.scanone$model <- rep("scanone", nrow(bayesint.result.scanone))

bayesint.result <- 
bayesint.result.cim %>% 
  full_join(bayesint.result.scanone, by = c("trait", "chr")) 

bayesint.result %>% dim() # 35 16 
```

```
## [1] 35 16
```

```r
for (i in 1:nrow(bayesint.result)){
  if(is.na(bayesint.result[i, "start.x"])){
    bayesint.result[i, c(3:9)] <- bayesint.result[i, c(10:16)]
  }
}

bayesint.result <- bayesint.result[,1:9] 
colnames(bayesint.result) <- gsub("\\.x$", "", colnames(bayesint.result)) 

bayesint.result %>% dim() # 35 9 
```

```
## [1] 35  9
```

```r
# make Table 2
length(bayesint.list.scanone) # 33
```

```
## [1] 33
```

```r
length(bayesint.list.cim) # 26 
```

```
## [1] 26
```

```r
bayesint.result.scanone <- as.tibble(bind_rows(bayesint.list.scanone,.id="trait")) %>% # combine list into tibble 
  dplyr::select(trait,chr,pos,markername,LOD) %>% 
  group_by(trait,chr) %>% 
  dplyr::summarize(start=min(pos, na.rm = T),end=max(pos, na.rm = T),pos = median(pos, na.rm = T), LOD=max(LOD))  

bayesint.result.scanone %>% dim() # 33 6 
```

```
## [1] 33  6
```

```r
bayesint.result.tmp <- as.tibble(bind_rows(bayesint.list.scanone,.id="trait")) %>% # combine list into tibble 
  dplyr::select(trait,chr,pos,markername,LOD) %>% 
  # separate(markername,into=c("chr1","Mbp"),sep="_", convert=TRUE) %>% 
  group_by(trait,chr) 

bayesint.result.tmp$index <- paste(bayesint.result.tmp$trait, bayesint.result.tmp$chr, bayesint.result.tmp$LOD)
bayesint.result.scanone$index <- paste(bayesint.result.scanone$trait, bayesint.result.scanone$chr, bayesint.result.scanone$LOD)

test <- 
bayesint.result.tmp %>% 
  anti_join(bayesint.result.scanone) %>% 
  dplyr::select(trait, chr, markername) %>% 
  mutate(index = paste(trait, chr, sep = "_")) 
```

```
## Joining, by = c("trait", "chr", "pos", "LOD", "index")
```

```r
tmp <- c()

for (i in seq_along(1:(nrow(test)/2))){
  tmp[i] <- paste(test$markername[i*2-1], test$markername[i*2-0], sep = "-") 
}

bayesint.result.scanone$index <- tmp

bayesint.result.scanone <- 
bayesint.result.scanone %>%
  mutate(start = round(start, 2), end = round(end, 2), pos = round(pos, 2), LOD = round(LOD, 2)) %>%
  unite(confidence_interval, start, end, sep = "-") %>% 
  mutate(flanking_marker = index) %>% 
  dplyr::select(-index)  

bayesint.result.scanone %>% dim() # 33 6 
```

```
## [1] 33  6
```

```r
### cim result 
bayesint.list.cim %>% length() # 26
```

```
## [1] 26
```

```r
bayesint.result <- as.tibble(bind_rows(bayesint.list,.id="trait")) %>% # combine list into tibble 
  dplyr::select(trait,chr,pos,markername,LOD) %>% 
  group_by(trait,chr) %>% 
  dplyr::summarize(start=min(pos, na.rm = T),end=max(pos, na.rm = T),pos = median(pos, na.rm = T), LOD=max(LOD)) 
  #for the high QTL peaks the interval width is 0.  That is overly precise and need to widen those.

bayesint.result %>% dim() # 26 6 
```

```
## [1] 26  6
```

```r
bayesint.result.tmp <- as.tibble(bind_rows(bayesint.list,.id="trait")) %>% # combine list into tibble 
  dplyr::select(trait,chr,pos,markername,LOD) %>% 
  # separate(markername,into=c("chr1","Mbp"),sep="_", convert=TRUE) %>% 
  group_by(trait,chr) 

bayesint.result.tmp$index <- paste(bayesint.result.tmp$trait, bayesint.result.tmp$chr, bayesint.result.tmp$LOD)
bayesint.result$index <- paste(bayesint.result$trait, bayesint.result$chr, bayesint.result$LOD)

test <- 
bayesint.result.tmp %>% 
  anti_join(bayesint.result) %>% 
  dplyr::select(trait, chr, markername) %>% 
  mutate(index = paste(trait, chr, sep = "_")) 
```

```
## Joining, by = c("trait", "chr", "pos", "LOD", "index")
```

```r
tmp <- c()

for (i in seq_along(1:(nrow(test)/2))){
  tmp[i] <- paste(test$markername[i*2-1], test$markername[i*2-0], sep = "-") 
}

bayesint.result$index <- tmp 

bayesint.result.cim <- 
bayesint.result %>% 
  mutate(start = round(start, 2), end = round(end, 2), pos = round(pos, 2), LOD = round(LOD, 2)) %>%
  unite(confidence_interval, start, end, sep = "-") %>% 
  mutate(flanking_marker = index) %>% 
  dplyr::select(-index)  

bayesint.result.cim %>% dim() # 26 6  
```

```
## [1] 26  6
```

```r
### combine cim & scanone result 
bayesint.result.scanone$model <- rep("scanone", nrow(bayesint.result.scanone))
bayesint.result.cim$model <- rep("cim", nrow(bayesint.result.cim))

bayesint.result.paper <- 
bayesint.result.cim %>% 
  full_join(bayesint.result.scanone, by = c("trait", "chr")) 

bayesint.result.paper %>% dim() # 32 12
```

```
## [1] 35 12
```

```r
colnames(bayesint.result.paper)
```

```
##  [1] "trait"                 "chr"                  
##  [3] "confidence_interval.x" "pos.x"                
##  [5] "LOD.x"                 "flanking_marker.x"    
##  [7] "model.x"               "confidence_interval.y"
##  [9] "pos.y"                 "LOD.y"                
## [11] "flanking_marker.y"     "model.y"
```

```r
for (i in 1:nrow(bayesint.result.paper)){
  if(is.na(bayesint.result.paper[i, "pos.x"])){
    bayesint.result.paper[i, c(3:7)] <- bayesint.result.paper[i, c(8:12)]
  }
}

bayesint.result.paper <- bayesint.result.paper[,1:7] 
colnames(bayesint.result.paper) <- gsub("\\.x$", "", colnames(bayesint.result.paper)) 

bayesint.result.paper %>% dim() # 35 7 
```

```
## [1] 35  7
```

```r
write.csv(bayesint.result.paper, file = "output/bayesint.result.paper.csv") 

### combine with allele effect information, allelic effect was calculated using fitqtl(), and values were added to bayesint.result.paper.csv  
bayesint.result.allele_effect <- read.csv("input/bayesint.result.paper_allele_effect.csv")

bayesint.result.paper.final <- 
bayesint.result.paper %>%
  left_join(bayesint.result.allele_effect, by = c("trait", "chr")) %>% 
  mutate(model = model.x) %>%
  dplyr::select(-model.x, -model.y, -X) 
```

```
## Warning: Column `trait` joining character vector and factor, coercing into
## character vector
```

```
## Warning: Column `chr` joining character vector and factor, coercing into
## character vector
```

```r
### add physical position interval  
bayesint.result.scanone <- as.tibble(bind_rows(bayesint.list.scanone,.id="trait")) %>% # combine list into tibble 
  dplyr::select(trait,chr,pos,markername,LOD) %>% 
  separate(markername,into=c("chr1","Mbp"),sep="_", convert=TRUE) %>% 
  group_by(trait,chr) %>% 
  dplyr::summarize(start=min(Mbp, na.rm = T),end=max(Mbp, na.rm = T),min_eQTL_LOD=min(LOD),max_eQTL_LOD=max(LOD))
```

```
## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 15 rows [5,
## 11, 26, 32, 35, 47, 53, 56, 62, 68, 83, 86, 92, 95, 98].
```

```r
bayesint.result.cim <- as.tibble(bind_rows(bayesint.list.cim,.id="trait")) %>% # combine list into tibble 
  dplyr::select(trait,chr,pos,markername,LOD) %>% 
  separate(markername,into=c("chr1","Mbp"),sep="_", convert=TRUE) %>% 
  group_by(trait,chr) %>% 
  dplyr::summarize(start=min(Mbp, na.rm = T),end=max(Mbp, na.rm = T),min_eQTL_LOD=min(LOD),max_eQTL_LOD=max(LOD))
```

```
## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 8 rows [2,
## 5, 29, 32, 44, 59, 68, 71].
```

```r
bayesint.result.cim$model <- rep("cim", nrow(bayesint.result.cim))
bayesint.result.scanone$model <- rep("scanone", nrow(bayesint.result.scanone))

bayesint.result <- 
bayesint.result.cim %>% 
  full_join(bayesint.result.scanone, by = c("trait", "chr")) 

for (i in 1:nrow(bayesint.result)){
  if(is.na(bayesint.result[i, "start.x"])){
    bayesint.result[i, c(3:7)] <- bayesint.result[i, c(8:12)]
  }
}

bayesint.result <- bayesint.result[,1:7] 
colnames(bayesint.result) <- gsub("\\.x$", "", colnames(bayesint.result)) 

bayesint.result %>% dim() # 35 7 
```

```
## [1] 35  7
```

```r
bayesint.result.physical <-
bayesint.result.paper.final %>% 
  left_join(bayesint.result, by= c("trait", "chr")) %>% 
  mutate(model = model.x) %>% 
  dplyr::select(-min_eQTL_LOD, -max_eQTL_LOD, -model.y, -model.x) 
  
write.csv(bayesint.result.physical, file = "output/bayesint.result.paper.physical.csv") 

# annotate QTL   
load("input/BnapusAnnotation.Rdata") 

traitQTL.annotated <- lapply(1:nrow(bayesint.result),function(row) { # for each trait/module 
  qtl <- bayesint.result[row,]  
  results <- subset(BnapusAnnotation, chrom==qtl$chr &
                    start >= qtl$start & # genes which fall into the QTL interval 
                    end <= qtl$end)
} 
)  

names(traitQTL.annotated) <- bayesint.result$trait 

traitQTL.annotated <- bind_rows(traitQTL.annotated,.id="trait") %>% # combine list into data.frame 
  mutate(chrom=as.character(chrom)) %>%
  left_join(bayesint.result,by=c("trait","chrom"="chr")) #get eQTL LOD

traitQTL.annotated <- 
traitQTL.annotated %>% 
  mutate(start = start.y, end = end.y) %>% 
  dplyr::select(-start.x, -end.x, -start.y, -end.y, -min_eQTL_LOD, -max_eQTL_LOD) 

traitQTL.annotated %>% dim() # 18647     11
```

```
## [1] 18647     9
```

```r
## get GO term for each gene
load("input/napus_GO_combined.Rdata") 

colnames(traitQTL.annotated)[3] <- "gene_ID"

traitQTL.annotated <- 
traitQTL.annotated %>% 
  left_join(napus_GO_combined) 
```

```
## Joining, by = "gene_ID"
```

```r
traitQTL.annotated %>% dim() #  18647    13 
```

```
## [1] 18647    11
```

```r
save(traitQTL.annotated, file =  "output/traitQTL.annotated.flipped.Rdata")  
```

### eQTL mapping result output  

```r
# determine cis- and trans- eQTL 


# output plot 


# trans-eQTL hotspot 
```

### integrate QTL and eQTL mapping result 

```r
# cis-coding candidates 


# cis-regulator candidates


# trans-regulator target candidates


# plot important candidates 
```

