# Li-eQTL-2018
Scripts for B. napus QTL and eQTL paper

This repository contains script, input, and output data for Li et al. 2018 on Brassica napus QTL and eQTL. 

The main scripts are the .MD and .RMD files in the scripts folder. 

* [phenotype-data-analysis.md](https://github.com/MaloofLab/Li-eQTL-2018/blob/master/scripts/phenotype-data-analysis.md) starts from raw phenotype data, plot histogram, check Pearson correlation between phenotypic traits, and calculate broad sense heritability

* [growth-model.Rmd](https://github.com/MaloofLab/Li-eQTL-2018/blob/master/scripts/growth-model.Rmd) includes the link to repository which does growth model fitting to the time-series data. The growth model uses R package brms to fit growth model, where the stats from current dataset were used as prior, multilevel model was applied to account for random effect from each F2 plant, and model was assessed using BIC/AIC, etc.   

* [gene-expression-analysis.md](https://github.com/MaloofLab/Li-eQTL-2018/blob/master/scripts/gene-expression-analysis.md) starts from RNAseq counts, TMM normalizes and fits glm model in edgeR for gene expression analysis. GO enrichment result of the differentially expressed genes is displayed in the file. 

* [SNP-calling.Rmd](https://github.com/MaloofLab/Li-eQTL-2018/blob/master/scripts/SNP-calling.Rmd) shows the SNP calling pipeline using Freebayes and GATK, SNP filtering of the output from both Feebayes and GATK, as well as how the two SNP datasets were combined. SNP annotation pipeline using snpEff is also included in this script. 

* [genetic-map-construction.Rmd](https://github.com/MaloofLab/Li-eQTL-2018/blob/master/scripts/genetic-map-construction.Rmd) starts from genotyping file of the F2 population, filters and formats the data, calculates two point recombination fraction, does clustering analysis to place markers on their respective linkage groups, and order markers on their linkage groups. 

* [genetic_map_construction_2.md](https://github.com/MaloofLab/Li-eQTL-2018/blob/master/scripts/genetic_map_construction_2.md) refines the genetic map by correcting double cross-over, using ripple function to correct linkage group A10, and removing problematic markers on the map. The final genetic map and the genetic VS. physical map comparison plots were displayed in this script. 

* [QTL-eQTL-mapping.Rmd](https://github.com/MaloofLab/Li-eQTL-2018/blob/master/scripts/QTL-eQTL-mapping.Rmd) does QTL mapping for all the phenotypic traits and eQTL mapping for all the expressed genes. 

* [QTL_eQTL_mapping_2.md](https://github.com/MaloofLab/Li-eQTL-2018/blob/master/scripts/QTL_eQTL_mapping_2.md) processes the QTL and eQTL mapping result to identify candidate genes and genetic regions associated with phenotypic traits. trans-eQTL hotspot analysis is also included in this script.      
