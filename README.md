# Li-eQTL-2018
Scripts for B. napus QTL and eQTL paper

This repository contains script, input, and output data for Li et al. 2018 on Brassica napus QTL and eQTL. 

The main scripts are the .MD and .RMD files in the scripts folder. 

* gene-expression-analysis.md starts from RNAseq counts, TMM normalizes and fits glm model in edgeR for gene expression analysis. GO enrichment result of the differentially expressed genes is displayed in the file. 

* SNP-calling.Rmd shows the SNP calling pipeline using Freebayes and GATK, SNP filtering of the output from both Feebayes and GATK, as well as how the two SNP datasets were combined. SNP annotation pipeline using snpEff is also included in this script. 

* genetic-map-construction.Rmd starts from genotyping file of the F2 population, filters and formats the data, calculates two point recombination fraction, does clustinerg analysis to place markers on their respective linkage groups, and order markers on their linkage groups. 

* genetic_map_construction_2.md refine the genetic map by correcting for double cross-over, using ripple function to correct linkage group A10, removing problematic markers on the map. The final genetic map and the genetic VS. physical map comparison plot were displayed in this script. 

* QTL-eQTL-mapping.Rmd does QTL mapping for all the phenotypic traits and eQTL mapping for all the expressed genes. 

* QTL_eQTL_mapping_2.md processes the QTL and eQTL mapping result to identify candidate genes and genetic regions for phenotypic traits. trans-eQTL hotspot analysis is also included in this script.      
