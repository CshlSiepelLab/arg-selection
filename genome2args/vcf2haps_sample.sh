#!/bin/bash
#$ -N vcf2haps_samp
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=4G

echo "_START_$(date)"
#module load Anaconda3/5.3.0

GITPATH='/grid/siepel/home_norepl/mo'
RELATE_PATH="/grid/siepel/home_norepl/mo/relate_v1.0.17_x86_64_static"

POP=$1 # e.g. CEU_99
GENES=$2

GENELIST=${GITPATH}/arg-selection/genome2args/${GENES}.txt
POPFILE=${POP}.txt
VCFDIR=${GITPATH}/${POP}_vcfs

while IFS=$'\t' read -r GENE CHR FROM TO; do
	mkdir $GENE
    echo Converting: ${POPFILE}_${GENE}_chr${CHR}_${FROM}_${TO}_gt.recode.vcf.gz > ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log
    ${RELATE_PATH}/bin/RelateFileFormats \
	--mode ConvertFromVcf \
	--haps ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}_dummy.haps \
	--sample ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.sample \
	-i ${VCFDIR}/${POPFILE}_${GENE}_chr${CHR}_${FROM}_${TO}_gt.recode \
	--chr $CHR >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log
    echo "RELATEFF_EXITSTAT_$?" >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log

	${GITPATH}/arg-selection/genome2args/vcf2hap.py ${VCFDIR}/${POPFILE}_${GENE}_chr${CHR}_${FROM}_${TO}_gt.recode.vcf.gz ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.haps >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log
	echo "VCF2HAP_EXITSTAT_$?" >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log

	${GITPATH}/arg-selection/genome2args/make_poplabels.py ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.sample CEU EUR ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}
	echo "MAKE_POPLABELS_EXITSTAT_$?" >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log

	${GITPATH}/arg-selection/genome2args/make_uniform_map.py ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.haps 1.25 ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}
	echo "MAKE_MAP_EXITSTAT_$?" >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log

done < $GENELIST

echo "_END_$(date)"
exit