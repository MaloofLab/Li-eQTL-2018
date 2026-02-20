# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains analysis scripts for Li et al. 2018, a *Brassica napus* QTL and eQTL paper. It is an R-based bioinformatics project organized as a sequential analysis pipeline.

## Running Scripts

All main analyses are R Markdown (`.Rmd`) files. Render them from the repo root:

```bash
Rscript -e "rmarkdown::render('scripts/<script-name>.Rmd')"
```

Plain `.R` scripts can be run directly:

```bash
Rscript scripts/<script-name>.R
```

**Important:** The `.Rmd` files use hardcoded `setwd("~/Desktop/F2_paper/submission/Li-eQTL-2018/scripts/")` calls. When running, either update this path or run from within the `scripts/` directory. All input paths are `../input/` and output paths are `../output/` relative to `scripts/`.

## Analysis Pipeline (Sequential Order)

The scripts must be run in this order, as each stage produces outputs consumed by later stages:

1. **`scripts/phenotype-data-analysis.Rmd`** — Processes raw growth trait data using growth model coefficients (from brms), calculates broad-sense heritability, Pearson correlations between phenotypic traits.

2. **`scripts/gene-expression-analysis.Rmd`** — Starts from `input/read.count.tsv`; TMM normalization via edgeR, VST transformation via DESeq2 (output: `output/vstMat.parent.Rdata`), GLM fitting, GO enrichment via GOseq.

3. **`scripts/SNP-calling.Rmd`** — SNP calling pipeline (Freebayes + GATK), filtering, merging, and snpEff annotation. The helper `scripts/SNP_info_extract_vcf_ann.sh` extracts mutation types from annotated VCF.

4. **`scripts/genetic-map-construction.Rmd`** → **`scripts/genetic_map_construction_2.Rmd`** — Builds genetic map from F2 genotyping data using onemap (two-point recombination fractions) and R/qtl (double crossover correction, ripple, marker filtering). Output: final genetic map CSV files in `output/`.

5. **`scripts/QTL-eQTL-mapping.Rmd`** — QTL mapping for all phenotypic traits and eQTL mapping for all expressed genes using R/qtl (`scanone`, `cim`), with permutation-based thresholds. Parallelized via snowfall. Large intermediate data stored on Box.

6. **`scripts/QTL_eQTL_mapping_2.Rmd`** — Post-processing of QTL/eQTL results: Bayesian confidence intervals, cis/trans classification, candidate gene identification.

7. **`scripts/trans_eQTL_hotspot_perm_transonly.R`** — Trans-eQTL hotspot permutation analysis using qtlhot.

8. **`scripts/scantwo/`** — Two-locus interaction scans (scantwo), parallelized across SLURM jobs.

## Key R Packages

- **qtl**: Core QTL/eQTL mapping (`scanone`, `cim`, `bayesint`, `read.cross`)
- **qtlhot**: Trans-eQTL hotspot analysis (`highlod`, `hotsize`)
- **snowfall**: Parallel computing for permutation analyses
- **edgeR** / **DESeq2**: RNA-seq normalization and differential expression
- **onemap**: Initial genetic map construction (two-point recombination)
- **GenomicRanges** / **GenomicFeatures** / **rtracklayer**: Genomic interval operations for cis/trans classification
- **goseq**: GO enrichment correcting for gene-length bias (uses `input/Bngo.list.Rdata`)
- **brms**: Bayesian growth model fitting (separate repository, outputs fed into `input/*.best.fitted.summary.csv`)

## Helper Utilities

`scripts/helper.R` provides shared functions sourced by other scripts:
- `GOseq.Bn.ORA()` — GO enrichment for *B. napus* gene lists
- `SNP.freebayes.reformat.Ae()` / `SNP.GATK.reformat()` — Parse VCF genotype fields
- `qtl_plot()` — ggplot2-based LOD score visualization faceted by chromosome
- `coef_extract()` / `coef_extract_width()` — Extract growth model parameters from brms output

## Data Organization

- `input/` — All raw and pre-processed input data (genotype files, phenotype CSVs, read counts, annotation, GO mappings)
- `output/` — Generated results (maps, QTL/eQTL results, normalized expression matrices, figures)
- Some large intermediate files (eQTL scanone results) were stored externally on Box and loaded via the `boxr` package

## Organism & Cross Design

- Species: *Brassica napus* (amphidiploid, A and C subgenomes, 19 chromosomes: A01–A10, C01–C09)
- Cross: F2 population from parents Da-Ae (Ae) × Da-Ol-1 (Ol)
- Genetic map: constructed from F2 genotyping data, markers ordered by MadMapper then refined in R/qtl
