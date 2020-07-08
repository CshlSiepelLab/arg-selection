#!/bin/bash
#$ -N genFeatures_array
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=8G

## Specify at submit time, match # of line in GENELIST file
# -t 1-13

echo "_START_$(date)"

module purge
module load Anaconda3/5.3.0

GITPATH='/sonas-hs/siepel/hpc_norepl/home/mo'
#GENELIST=${GITPATH}/arg-selection/genome2args/pos_sel_genes.txt
GENELIST=$1
MINDAF=$2

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
${GITPATH}/arg-selection/genome2args/genFeatures.py ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}_wg.trees ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}_DAF${MINDAF} 2 $MINDAF
echo "_EXITSTAT_$?"

echo "_END_$(date)"
exit