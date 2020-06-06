#!/bin/bash
#$ -N pkl2fea
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=16G

## These should be passed in while submitting the job
# -t 1-60
# -tc 30

echo "_START_$(date)"

module purge
module load Anaconda3/5.3.0

GITPATH='/sonas-hs/siepel/hpc_norepl/home/mo'

ZEROBASE=$1 # Flag indicating if threads are zero-based (1- zero base; 0- NOT zero base)
HANDLE=$2
TAG=$3 # some unique identifier

#FLANKINGTR=2 <== this is part of the .py code now

${GITPATH}/arg-selection/sim2args/pkl2fea.py $HANDLE $((SGE_TASK_ID-ZEROBASE)) $TAG
echo "_EXITSTAT_$?"

echo "_END_$(date)"
exit
