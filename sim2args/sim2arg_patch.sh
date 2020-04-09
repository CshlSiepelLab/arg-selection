#!/bin/bash
#$ -N sim2arg_patch
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=48G

# -pe threads 32

echo "_START_$(date)"

module purge
module load Anaconda3/5.3.0

GITPATH='/sonas-hs/siepel/hpc_norepl/home/mo'

ZEROBASE=$1 # Flag indicating if threads are zero-based (1- zero base; 0- NOT zero base)
HANDLE=$2
TAG=$3 # some unique identifier
THREAD=$4

NE=188088
FLANKINGTR=5

echo Patching ${HANDLE}_${TAG}_${THREAD}
${GITPATH}/arg-selection/sim2args/sim2arg.py $FLANKINGTR $HANDLE $NE $((THREAD-ZEROBASE)) $TAG
echo "_EXITSTAT_$?"

echo "_END_$(date)"
exit
