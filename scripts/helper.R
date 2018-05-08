# 1) GO annotaion 

# ORA with GOseq (Brassica napa version)
library(ShortRead);library(goseq);library(GO.db);library("annotate")
library(WGCNA);library(ggplot2);library(reshape2);library(scales); library (plyr)

Bn_cdna<-readDNAStringSet("input/Brassica_napus.annotation_v5.cds_modified.fa")
bias<-nchar(Bn_cdna)
names(bias)<-names(Bn_cdna)

#  bias.data vector must have the same length as DEgenes vector!
# convert to list (onetime procedure)
# Bngo<-read.table("input/Brassica_napus_GO",header=FALSE) 
# 
# Bngo.list <- tapply(as.character(Bngo$V2),Bngo$V1,c)
# save(Bngo.list,file="input/Bngo.list.Rdata") 
# # Using manually entered categories.
# # Calculating the p-values...
# # 'select()' returned 1:1 mapping between keys and columns
# Bngo.DF<-as.data.frame(Bngo.list)
# Bngo.DF$gene<-rownames(Bngo.DF)
# do.call(rbind.data.frame, Bngo.list)
# Bngo.DF2<-do.call(rbind.data.frame,Bngo.list) 
# Bngo.DF3 <- ldply (Bngo.list, data.frame)
# names(Bngo.DF3)<-c("gene","GO") #if Brgo.list does not work in goseq, use DF3.

load("input/Bngo.list.Rdata")

GOseq.Bn.ORA<-function(genelist,padjust=0.05,ontology="BP") { # return GO enrichment table, padjus, padjust=0.05
  TF<-(names(bias) %in% genelist)*1
  names(TF)<-names(bias) 
  pwf<-nullp(TF,bias.data=bias) 
  GO.pval <- goseq(pwf,gene2cat=Bngo.list,use_genes_without_cat=TRUE) 
  
  if(ontology=="BP") {
    GO.pval2<-subset(GO.pval,ontology=="BP")
  } else if(ontology=="CC") {
    GO.pval2<-subset(GO.pval,ontology=="CC")
  } else {
    GO.pval2<-subset(GO.pval,ontology=="MF")
  }
  
  GO.pval2$over_represented_padjust<-p.adjust(GO.pval2$over_represented_pvalue,method="BH")
  if(GO.pval2$over_represented_padjust[1]>padjust) stop("no enriched GO")
  else {
    enriched.GO<-GO.pval2[GO.pval2$over_represented_padjust<padjust,]
    print("enriched.GO is")
    print(enriched.GO)
    
    ## write Term and Definition
    for(i in 1:dim(enriched.GO)[1]) {
      enriched.GO$Term[i]<-Term(GOTERM[[enriched.GO[i,"category"]]])
      enriched.GO$Definition[i]<-Definition(GOTERM[[enriched.GO[i,"category"]]])
    }
    return(enriched.GO)
  }
}  
