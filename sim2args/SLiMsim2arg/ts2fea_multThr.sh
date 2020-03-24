#!/bin/bash
#$ -N ts2Feature_arr
#$ -S /bin/bash
#$ -cwd
#$ -t 1-200
#$ -tc 200
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=8G

echo "_START_$(date)"

module purge
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6

GITPATH="/sonas-hs/siepel/hpc_norepl/home/mo/arg-selection"

MODE=$1 # <`n`/`s`>
META=$2 # <no_sims/meta_file_path>

THR=$SGE_TASK_ID
TOTTHR=200

INPREF=$3
TTYPE=$4 # <`tru`/`inf`>
OUTPREF=$5

${GITPATH}/sim2args/SLiMsim2arg/ts2feature.py $MODE $META $THR $TOTTHR $INPREF $TTYPE $OUTPREF
echo "_EXITSTAT_$?"

echo "_END_$(date)"
exit
