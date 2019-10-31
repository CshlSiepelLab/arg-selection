#!/bin/bash

echo "_START_$(date)"

module purge
module load GCC/7.3.0-2.30

POPFILE=$1
CHR=$2 # note that non autosomes have different file naming from autosomes
FROM=$3
TO=$4
GENE=$5

DATAPATH=/usr/data/1000G
OUTPATH=/sonas-hs/siepel/hpc_norepl/home/mo/sel_coef_empirical/vcfs
# --from-bp <integer> 
# --to-bp <integer>
# These options specify a lower bound and upper bound for a range of sites to be processed. Sites with positions less than or greater than these values will be excluded.

# --keep <filename>
# Provide files containing a list of individuals to either include or exclude in subsequent analysis. Each individual ID (as defined in the VCF headerline) should be included on a separate line. No header line is expected.

vcftools --vcf ${DATAPATH}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf --keep $POPFILE --remove-indels --min-alleles 2 --max-alleles 2 --chr $CHR --from-bp $FROM --to-bp $TO --recode --recode-INFO AA --out ${OUTPATH}/${POPFILE}_${GENE}_chr${CHR}_${FROM}_${TO}_gt

bgzip ${OUTPATH}/${POPFILE}_${GENE}_chr${CHR}_${FROM}_${TO}_gt.recode.vcf

tabix -p vcf ${OUTPATH}/${POPFILE}_${GENE}_chr${CHR}_${FROM}_${TO}_gt.recode.vcf.gz

echo "_EXITSTAT_$?"
echo "_END_$(date)"
