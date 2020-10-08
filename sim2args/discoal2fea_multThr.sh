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

#module purge
#module load Anaconda3/5.3.0

GITPATH='/grid/siepel/home_norepl/mo'
EXEC='discoal2fea.py' # inferred pipeline
#EXEC='discoal2fea_neu_bugfix.py' # fixing bug in neutral sim site selection
#EXEC='discoalT2fea.py' # true trees from discoal

HANDLE=$1 # e.g. "soft_lu_110_swp"
NO_PKL=$2
NO_THR=$SGE_TASK_LAST

${GITPATH}/arg-selection/sim2args/${EXEC} $HANDLE $SGE_TASK_ID $NO_PKL $NO_THR
echo "_EXITSTAT_$?"

echo "_END_$(date)"
exit