#!/bin/bash

# goal: get SNP mutation type from annotated vcf file 

SNP=$1
vcf=$2

# file_in format
#    CHROM      POS REF ALT Ae.gt Ol.gt
# 1 chrA08 10022021   G   C   1/1   0/0
# 2 chrA08 10030368   A   G   1/1   0/0
# 3 chrA08 10030746   T   G   1/1   0/0
# 4 chrA08 10030784   C   A   1/1   0/0
# 5 chrA08 10034337   T   A   1/1   0/0
# 6 chrA08 10034794   G   T   1/1   0/0

nl=`cat $SNP | wc -l`
nl=`expr $nl - 1` 
cat $SNP | sed 's/"//g' | sed 's/,/ /g' | awk '{print $2 "_" $3}' | sed 's/chr/S/g' | tail -$nl  | sed 's/S//g' | awk '{print "chr" $0}' > significant_SNP

awk 'NR==FNR{a[$1]=$0}NR>FNR{if ($1 in a) print $0}' significant_SNP $vcf  | awk '{print $2 "|" $3 "|" $5 "|" $6 "|" $9 "|" $10}' | awk 'BEGIN{FS="|"}{print $1, $2, $3, $4, $6, $7}' > significant_SNP_annotation

rm significant_SNP
 
