#!/bin/bash
#$ -N genFeatures_array
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=10G

## Specify at submit time, match # of line in GENELIST file
# -t 1-16

echo "_START_$(date)"

# module purge
# module load Anaconda3/5.3.0

GITPATH='/grid/siepel/home_norepl/mo'
GENES=$1
MINDAF=$2

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

# usage: $./genFeatures.py <.trees PATH> <out prefix> <no_ft> <min_DAF>

echo Feature Extraction ${GENE}_chr${CHR}_${FROM}_${TO}
${GITPATH}/arg-selection/genome2args/genFeatures.py ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}_wg.trees ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.haps.meta ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}_DAF${MINDAF}_gt 2 $MINDAF
echo "_EXITSTAT_$?"

echo "_END_$(date)"
exit
