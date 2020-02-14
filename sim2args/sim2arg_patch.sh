#!/bin/bash
#$ -N sim2arg_patch
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=32G

# -pe threads 32

echo "_START_$(date)"

module purge
module load Anaconda3/5.3.0

GITPATH='/sonas-hs/siepel/hpc_norepl/home/mo'

HANDLE=$1
TAG=$2 # some unique identifier
THREAD=$3

NE=188088
FLANKINGTR=2

echo Patching ${HANDLE}_${TAG}_${THREAD}
${GITPATH}/arg-selection/sim2args/sim2arg.py $FLANKINGTR $HANDLE $NE $THREAD $TAG
echo "_EXITSTAT_$?"

echo "_END_$(date)"
exit
