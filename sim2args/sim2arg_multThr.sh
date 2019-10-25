#!/bin/bash
#$ -N sim2arg
#$ -S /bin/bash
#$ -cwd
#$ -t 1-50
#$ -tc 25
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=48G

echo "_START_$(date)"

module purge
module load Anaconda3/5.3.0

#MODE=$1
HANDLE=$1
NE=$2
#TAG=$4 # some unique identifier

#        usage: $./sim2arg.py <fea> <handle> <Ne> <thread #, 1-based> <TAG>
#            Takes a *partitioned* pickle file, runs RELATE to infer ARGs and extract features.
#            <fea> determines the features extracted:
#                1: single gene tree at site of interest
#                3: incorporate the immediate surrounding trees
#                5: incorporate the two immediate surrounding trees

#./sim2arg.py $MODE discoal_${HANDLE} $NE $SGE_TASK_ID ${MODE}_${TAG}
./sim2arg_neu.py 1 discoal_${HANDLE} $NE $SGE_TASK_ID OCT2
#./sim2arg_RENTplus.py discoal_${HANDLE} $SGE_TASK_ID $TAG

echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit
