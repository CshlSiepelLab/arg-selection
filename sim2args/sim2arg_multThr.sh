#!/bin/bash
#$ -N sim2arg
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=32G

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
#RECOMB=$4 # recombination rate (in 1e-8 unit)

NE=188088
FLANKINGTR=2

# usage: $./sim2arg.py <no_st> <handle> <Ne> <thread #, 1-based> <TAG>
#     Takes a *partitioned* pickle file, runs RELATE to infer ARGs and extract features.
#     - <no_st> >=0 : number of flanking gene trees to include on EACH side for feature extraction

${GITPATH}/arg-selection/sim2args/sim2arg.py $FLANKINGTR $HANDLE $NE $((SGE_TASK_ID-ZEROBASE)) $TAG
echo "_EXITSTAT_$?"

echo "_END_$(date)"
exit
