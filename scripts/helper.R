# 1) GO annotaion 

# ORA with GOseq (Brassica napa version)
library(ShortRead);library(goseq);library(GO.db);library("annotate")
library(WGCNA);library(ggplot2);library(reshape2);library(scales); library (plyr)

Bn_cdna<-readDNAStringSet("../input/Brassica_napus.annotation_v5.cds_modified.fa")
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

load("../input/Bngo.list.Rdata")

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

SNP.freebayes.reformat.Ae <- function(vcf, vcf.header){ 
  colnames(vcf) <- vcf.header
  head(vcf)
  
  vcf$Ae[is.na(vcf$Ae)] <- "NA:NA:NA:NA:NA:NA:NA"
  
  Ae.tmp.unique <- matrix(
    unlist(strsplit(vcf$Ae,split = ":")),
    nrow=nrow(vcf),  
    byrow=TRUE
  ) 
  
  colnames(Ae.tmp.unique) <- paste("Ae",c("gt","tot.depth","ref.depth","ref.qual","alt.depth","alt.qual","gen.lik"),sep="_")
  
  vcf.reform <- cbind(vcf,Ae.tmp.unique,stringsAsFactors=FALSE)
  
  vcf.reform[,c("Ae_tot.depth","Ae_ref.depth","Ae_ref.qual","Ae_alt.depth","Ae_alt.qual")] <- 
    apply(vcf.reform[,c("Ae_tot.depth","Ae_ref.depth","Ae_ref.qual","Ae_alt.depth","Ae_alt.qual")],
          2,
          as.numeric) 
  
  return(vcf.reform)  
} 

# function 
SNP.GATK.reformat <- function(vcf, vcf.header){ # input are vcf file from GATK after filtering & header of the vcf file 
  colnames(vcf) <- vcf.header 
  
  # correct the file 
  problem.row <- which(vcf$FORMAT!="GT:AD:DP:GQ:PL")
  vcf.corrected <- vcf[-c(problem.row),]
  
  # Before splitting add NAs to blank cells 
  vcf.corrected$Ae[vcf.corrected$Ae=="./.:.:.:.:."] <- "NA:NA,NA:NA:NA:NA,NA,NA"
  Ae.GATK <- matrix(
    unlist(strsplit(vcf.corrected$Ae,split = ":")), 
    nrow=nrow(vcf.corrected),  
    byrow=TRUE
  ) 
  colnames(Ae.GATK) <- paste("Ae",c("gt","ref.alt.depth","approx.depth","genotype.qual","Phred.score"),sep="_")
  
  vcf.corrected$Ol[vcf.corrected$Ol=="./.:.:.:.:."] <- "NA:NA,NA:NA:NA:NA,NA,NA" 
  Ol.GATK <- matrix(
    unlist(strsplit(vcf.corrected$Ol,split = ":")),
    nrow=nrow(vcf.corrected),  
    byrow=TRUE
  ) 
  colnames(Ol.GATK) <- paste("Ol",c("gt","ref.alt.depth","approx.depth","genotype.qual","Phred.score"),sep="_")
  vcf.reform <- cbind(vcf.corrected, Ae.GATK, Ol.GATK,stringsAsFactors=FALSE) 
  
  vcf.reform[,c("Ae_approx.depth","Ae_genotype.qual",
                "Ol_approx.depth","Ol_genotype.qual")] <-
    apply(vcf.reform[,c("Ae_approx.depth","Ae_genotype.qual",
                        "Ol_approx.depth","Ol_genotype.qual")],
          2,
          as.numeric
    )
  return(vcf.reform)  
}   
