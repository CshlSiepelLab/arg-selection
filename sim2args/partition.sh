#!/bin/bash
#$ -N part_pkl
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=16G

## Specify at submit time
# -t 1-10
# -tc 10

echo "_START_$(date)"

module purge
module load Anaconda3/5.3.0

GITPATH='/sonas-hs/siepel/hpc_norepl/home/mo'
FILEPATH=$1 #path to pickle files (entire prefix, w/o trailing '_')
PKLPREF=$2

PARTITIONI=$((SGE_TASK_ID-1))
PARTOPI=$3

${GITPATH}/arg-selection/sim2args/partition.py $FILEPATH $PKLPREF $PARTITIONI $PARTOPI
echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit
