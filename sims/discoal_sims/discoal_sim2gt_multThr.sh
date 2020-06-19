#!/bin/bash
#$ -N discoal_sim
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=8G

## specify at submit time
# -t 1-3
# -tc 3

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.7.4

GITPATH='/sonas-hs/siepel/hpc_norepl/home/mo'

echo "_START_$(date)"

COND=$1 # ("swp" or "neu")
HDL=$2
SC_RANGE=$3
START_FREQ_RANGE=$4
END_FREQ_RANGE=$5
NOCHR=$6
NOSIMS_PTHR=$7

${GITPATH}/arg-selection/sims/discoal_sims/discoal_sim2gt.py $COND $HDL $SC_RANGE $START_FREQ_RANGE $END_FREQ_RANGE $NOCHR $SGE_TASK_ID $NOSIMS_PTHR
echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit
