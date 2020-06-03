#!/bin/bash
#$ -N discoal2fea
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=8G

## These should be passed in while submitting the job
# -t 1-60
# -tc 30

echo "_START_$(date)"

module purge
module load Anaconda3/5.3.0

GITPATH='/sonas-hs/siepel/hpc_norepl/home/mo'

MODE=$1 # "-n" or "-s"
HANDLE=$2 # e.g. "pw_32_swp"
NOREPL=$3
NOPART=$4

${GITPATH}/arg-selection/sim2args/discoal2fea.py $MODE $HANDLE $NOREPL $NOPART $((SGE_TASK_ID-1))
echo "_EXITSTAT_$?"

echo "_END_$(date)"
exit
