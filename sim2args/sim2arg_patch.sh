#!/bin/bash
#$ -N sim2arg
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=48G


echo "_START_$(date)"

module purge
module load Anaconda3/5.3.0

MODE=$1
HANDLE=$2
NE=$3
THREAD=$4
TAG=$5

echo Patching ${HANDLE}_${MODE}_${TAG}_${THREAD}
./sim2arg.py $MODE discoal_${HANDLE} $NE $THREAD ${MODE}_${TAG}

echo "_EXITSTAT_$?"
echo "_END_$(date)"

exit
