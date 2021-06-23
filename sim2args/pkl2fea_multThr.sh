#!/bin/bash
#$ -N pkl2fea
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=8G

## These should be passed in while submitting the job
# -t 1-60
# -tc 30

echo "_START_$(date)"

# module purge
# module load EBModules
# module load Anaconda3/2019.10
# module load R/3.6.2-foss-2019b
# DO THIS FIRST: $ conda activate local

GITPATH='/grid/siepel/home_norepl/mo'

# ZEROBASE=$1 # Flag indicating if threads are zero-based (1- zero base; 0- NOT zero base)
# HANDLE=$2
# TAG=$3 # some unique identifier
# #FLANKINGTR=2 # <== this is part of the .py code now
# ${GITPATH}/arg-selection/sim2args/pkl2fea.py $HANDLE $((SGE_TASK_ID-ZEROBASE)) $TAG

PKL_PREF=$1
TOT_FILES=$2
OUT_HANDLE=$3 # some unique identifier
${GITPATH}/arg-selection/sim2args/pkl2fea_noprep.py $PKL_PREF $((SGE_TASK_ID-1)) $SGE_TASK_LAST $TOT_FILES $OUT_HANDLE

echo "_EXITSTAT_$?"

echo "_END_$(date)"
exit
