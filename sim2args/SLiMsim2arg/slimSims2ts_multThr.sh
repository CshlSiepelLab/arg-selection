#!/bin/bash
#$ -N slim2ts_arr
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=32G

## Specify with qsub
# -t 1-50
# -tc 50

echo "_START_$(date)"

module purge
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load GSL/2.5

GITPATH="/sonas-hs/siepel/hpc_norepl/home/mo/arg-selection"

PARAMF=$1
SCALE=$2
#MODE=$3 # <`n`/`s`>
#META=$4 # <no_sims/meta_file_path>
NOSIMS=$3
THR=$SGE_TASK_ID
TOTTHR=$4
INPREF=$5
OUTPREF=$6

${GITPATH}/sim2args/SLiMsim2arg/slimSims2ts.py $PARAMF $SCALE $NOSIMS $THR $TOTTHR $INPREF $OUTPREF
echo "_EXITSTAT_$?"

echo "_END_$(date)"
exit
