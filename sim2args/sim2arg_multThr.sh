#!/bin/bash
#$ -N sim2arg
#$ -S /bin/bash
#$ -cwd
#$ -t 1-50
#$ -tc 50
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=40G

echo "_START_$(date)"

module purge
module load Anaconda3/5.3.0

GITPATH='/sonas-hs/siepel/hpc_norepl/home/mo'

HANDLE=$1
TAG=$2 # some unique identifier

NE=188088
FLANKINGTR=2

# usage: $./sim2arg.py <no_st> <handle> <Ne> <thread #, 1-based> <TAG>
#     Takes a *partitioned* pickle file, runs RELATE to infer ARGs and extract features.
#     - <no_st> >=0 : number of flanking gene trees to include on EACH side for feature extraction

${GITPATH}/arg-selection/sim2args/sim2arg.py $FLANKINGTR $HANDLE $NE $SGE_TASK_ID $TAG
echo "_EXITSTAT_$?"

echo "_END_$(date)"
exit
