#!/bin/bash
#$ -N genFeatures_array
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=10G

## Specify at submit time
# -t 1-16

echo "_START_$(date)"

# module purge
# module load Anaconda3/5.3.0

GITPATH='/grid/siepel/home_norepl/mo'

INTAG=$1
MINMAXDAF=$2
OUTPREF=$3
CHR=$4

THR=$((SGE_TASK_ID-1))
NO_THR=$SGE_TASK_LAST

echo Feature Extraction ${INTAG}_chr${CHR} Thread $THR
${GITPATH}/arg-selection/genome2args/genFeatures_win.py ${INTAG}_chr${CHR}/${INTAG}_${CHR}_wg.trees ${INTAG}_chr${CHR}/${INTAG}_${CHR}.meta ${OUTPREF}_${INTAG}_chr${CHR} $MINMAXDAF $NO_THR $THR
echo "_EXITSTAT_$?"

echo "_END_$(date)"
exit