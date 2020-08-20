#!/bin/bash
#$ -N emp_pipeline
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=16G

echo "_START_$(date)"

GITPATH='/grid/siepel/home_norepl/mo'
RELATE_PATH="/grid/siepel/home_norepl/mo/relate_v1.0.17_x86_64_static"
#GENELIST=${GITPATH}/arg-selection/genome2args/pos_sel_genes.txt

GENE=$1
CHR=$2
FROM=$3
TO=$4

cd $GENE

### ARG INFERENCE ###
echo Running Relate on ${GENE}_chr${CHR}_${FROM}_${TO}

${RELATE_PATH}/bin/Relate \
	--mode All \
	-m 2.5e-8 \
	-N 376176 \
	--haps ${GENE}_chr${CHR}_${FROM}_${TO}.haps \
	--sample ${GENE}_chr${CHR}_${FROM}_${TO}.sample \
	--map ${GENE}_chr${CHR}_${FROM}_${TO}.map \
	-o ${GENE}_chr${CHR}_${FROM}_${TO}
echo "RELATE_EXITSTAT_$?"

# re-estimate branch lengths
${RELATE_PATH}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
            -i ${GENE}_chr${CHR}_${FROM}_${TO} \
            -m 2.5e-8 \
            --poplabels ${GENE}_chr${CHR}_${FROM}_${TO}.poplabels \
            --threshold 10 \
            -o ${GENE}_chr${CHR}_${FROM}_${TO}_popsize
echo "POPSIZE_EXITSTAT_$?"

# re-estimate branch length for ENTIRE genealogy
${RELATE_PATH}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
            -i ${GENE}_chr${CHR}_${FROM}_${TO} \
            -m 2.5e-8 \
            --poplabels ${GENE}_chr${CHR}_${FROM}_${TO}.poplabels \
            --threshold 0 \
            --coal ${GENE}_chr${CHR}_${FROM}_${TO}_popsize.coal \
            --num_iter 1 \
            -o ${GENE}_chr${CHR}_${FROM}_${TO}_wg
echo "WG_EXITSTAT_$?"

${RELATE_PATH}/bin/RelateFileFormats \
	--mode ConvertToTreeSequence \
	-i ${GENE}_chr${CHR}_${FROM}_${TO}_wg \
	-o ${GENE}_chr${CHR}_${FROM}_${TO}_wg
echo "TS_CONVERSION_EXITSTAT_$?"

### FEATURE EXTRACTION ###
MINDAF=0.05

echo Feature Extraction ${GENE}_chr${CHR}_${FROM}_${TO}
${GITPATH}/arg-selection/genome2args/genFeatures.py ${GENE}_chr${CHR}_${FROM}_${TO}_wg.trees ${GENE}_chr${CHR}_${FROM}_${TO}_DAF${MINDAF} 2 $MINDAF
echo "FEA_EXITSTAT_$?"

echo "_END_$(date)"
exit
