#!/bin/bash
#$ -N RELATE_array
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=8G

## Specify at submit time, match # of line in GENELIST file
# -t 1-16

echo "_START_$(date)"

GITPATH='/grid/siepel/home_norepl/mo'
RELATE_PATH="/grid/siepel/home_norepl/mo/relate_v1.0.17_x86_64_static"

GENES=$1
GENELIST=${GITPATH}/arg-selection/genome2args/${GENES}.txt

GENEL=($(awk '{print $1}' $GENELIST))
CHRL=($(awk '{print $2}' $GENELIST))
FROML=($(awk '{print $3}' $GENELIST))
TOL=($(awk '{print $4}' $GENELIST))

#for ((SGE_TASK_ID=1;SGE_TASK_ID<=13;SGE_TASK_ID++)); do
IDX=$((SGE_TASK_ID-1))
GENE=${GENEL[$IDX]}
CHR=${CHRL[$IDX]}
FROM=${FROML[$IDX]}
TO=${TOL[$IDX]}

cd $GENE
echo Running Relate on ${GENE}_chr${CHR}_${FROM}_${TO}

#../relate_v1.0.16_x86_64_static/scripts/PrepareInputFiles/PrepareInputFiles.sh \
#	--haps ${GENE}_chr${CHR}_${FROM}_${TO}.haps \
#	--sample ${GENE}_chr${CHR}_${FROM}_${TO}.sample \
#	--ancestor ../ancestral_fasta/homo_sapiens_ancestor_${CHR}.fa \
#	-o ${GENE}_chr${CHR}_${FROM}_${TO}_input

#gunzip ${GENE}_chr${CHR}_${FROM}_${TO}_input.haps.gz
#gunzip ${GENE}_chr${CHR}_${FROM}_${TO}_input.sample.gz

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

echo "_END_$(date)"
exit
