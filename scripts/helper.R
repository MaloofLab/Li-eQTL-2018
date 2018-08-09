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

# function to plot qtl result 
qtl_plot <- function(input,              # data frame input from scanone
                     mult.pheno = FALSE, # multiple phenotypes?
                     model = "normal",   # model used in scanone
                     chrs = NA,          # chromosomes to display
                     lod = NA,           # LOD threshold
                     rug = FALSE,        # plot marker positions as rug?
                     ncol = NA,          # number of columns for facetting
                     labels = NA,         # optional dataframe to plot QTL labels
                     title = NA          # title of the figure 
) {
  
  # if we have multiple phenotypes and/or a 2part model, gather input
  if (mult.pheno & model == "2part") {
    input <- gather(input, group, lod, grep("pheno", colnames(input)))
  } else if (mult.pheno) {
    input <- gather(input, group, lod, grep("pheno", colnames(input)))
  } else if (model == "2part") {
    input <- gather(input, method, lod, lod.p.mu:lod.mu)
  }
  
  # if not all chromosomes should be displayed, subset input
  if (!is.na(chrs)[1]) {
    input <- input[as.character(input$chr) %in% chrs, ]
  }
  
  # if there is more than one LOD column, gather input
  if (!any(colnames(input) == "lod")) {
    input$lod <- input[, grep("lod", colnames(input))]
  }
  
  # if no number of columns for facetting is defined, plot all in one row
  if (is.na(ncol)) {
    ncol <- length(unique(input$chr))
  }
  
  # if labels are set and there is no name column, set from rownames
  if (!is.na(labels)[1]) {
    if (is.null(labels$name)) {
      labels$name <- rownames(labels)
    }
  }
  
  # if tile are given
  if (!is.na(title)){
    labs(title = "title")
  }
  
  # plot input data frame position and LOD score
  plot <- ggplot(input, aes(x = pos, y = lod)) + {
    
    # if LOD threshold is given, plot as horizontal line
    if (!is.na(lod)[1] & length(lod) == 1) geom_hline(yintercept = lod, linetype = "dashed")
  } + {
    
    if (!is.na(lod)[1] & length(lod) > 1) geom_hline(data = lod, aes(yintercept = lod, linetype = group))
  } + {
    
    # plot rug on bottom, if TRUE
    if (rug) geom_rug(size = 0.1, sides = "b")
  } + {
    
    # if input has column method but not group, plot line and color by method
    if (!is.null(input$method) & is.null(input$group)) geom_line(aes(color = method), size = 1, alpha = 0.6)
  } + {
    
    # if input has column group but not method, plot line and color by group
    if (!is.null(input$group) & is.null(input$method)) geom_line(aes(color = group), size = 1, alpha = 0.6)
  } + { 
    
    # if input has columns method and group, plot line and color by method & linetype by group
    if (!is.null(input$group) & !is.null(input$method)) geom_line(aes(color = method, linetype = group), size = 1, alpha = 0.6)
  } + {
    
    # set linetype, if input has columns method and group
    if (!is.null(input$group) & !is.null(input$method)) scale_linetype_manual(values = c("solid", "twodash", "dotted"))
  } + { 
    # if input has neither columns method nor group, plot black line
    if (is.null(input$group) & is.null(input$method)) geom_line(size = 1, alpha = 0.6)
  } + {
    
    # if QTL positions are given in labels df, plot as point...
    if (!is.na(labels)[1]) geom_point(data = labels, aes(x = pos, y = lod))
  } + {
    
    # ... and plot name as text with ggrepel to avoid overlapping
    if (!is.na(labels)[1]) geom_text_repel(data = labels, aes(x = pos, y = lod, label = name), nudge_y = 0.5)
  } + 
    # facet by chromosome
    facet_wrap(~ chr, ncol = ncol, scales = "free_x") +
    # minimal plotting theme
    theme(panel.background = element_blank()) +
    # increase strip title size
    theme(strip.text = element_text(face = "bold", size = 8)) +
    # use RcolorBrewer palette
    scale_color_brewer(palette = "Set1") + ### default Set1 gave the best color combination  
    # scale_color_brewer(palette = "Set3") + ### 
    # no plot 
    theme(legend.position = "bottom") + ### this can be changed or , axis.text.x=element_blank() for making supplementary figures
    # set font size to 8 
    theme(text = element_text(size=8)) + 
    theme(axis.text=element_text(size=rel(0.8))) + 
    # Change plot labels
    labs(x = "Chromosome", 
         y = "LOD",
         color = "",
         linetype = "", 
         title = title)
  
  print(plot)
} 

coef_extract <- function(model_summary, ID){
  coef_Hmax <- model_summary[grep("r_ID__Hmax", model_summary$X) ,"mean"]
  coef_k <- model_summary[grep("r_ID__k", model_summary$X) ,"mean"] 
  coef_I <- model_summary[grep("r_ID__I", model_summary$X) ,"mean"] 
  intercept_Hmax <- model_summary[model_summary$X == "b_Hmax_Intercept", "mean"]
  intercept_k <- model_summary[model_summary$X == "b_k_Intercept", "mean"]
  intercept_I <- model_summary[model_summary$X == "b_I_Intercept", "mean"]
  
  growth_model_trait <- data.frame(line_ID = ID,
                                   Hmax = coef_Hmax + intercept_Hmax,
                                   k = coef_k + intercept_k,
                                   I = coef_I + intercept_I
  )
  return(growth_model_trait)
}

coef_extract_width <- function(model_summary, ID){
  
  model_summary <- width
  
  coef_Hmax <- model_summary[grep("r_ID__Hmax", model_summary$X) ,"mean"]
  coef_k <- model_summary[grep("r_ID__k", model_summary$X) ,"mean"] 
  coef_delta <- model_summary[grep("r_ID__delta", model_summary$X) ,"mean"] 
  intercept_Hmax <- model_summary[model_summary$X == "b_Hmax_Intercept", "mean"]
  intercept_k <- model_summary[model_summary$X == "b_k_Intercept", "mean"]
  intercept_delta <- model_summary[model_summary$X == "b_delta_Intercept", "mean"]
  
  growth_model_trait <- data.frame(line_ID = ID,
                                   Hmax = coef_Hmax + intercept_Hmax,
                                   k = coef_k + intercept_k,
                                   delta = coef_delta + intercept_delta
  )
  return(growth_model_trait)
} 