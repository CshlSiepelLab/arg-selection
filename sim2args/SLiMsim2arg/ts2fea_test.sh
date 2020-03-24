#!/bin/bash
#$ -N ts2Feature
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=32G

echo "_START_$(date)"

module purge
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6

GITPATH="/sonas-hs/siepel/hpc_norepl/home/mo/arg-selection"

MODE=$1 # <`n`/`s`>
META=$2 # <no_sims/meta_file_path>
THR=$3
TOTTHR=$4
INPREF=$5
TTYPE=$6 # <`tru`/`inf`>
OUTPREF=$7

${GITPATH}/sim2args/SLiMsim2arg/ts2feature.py $MODE $META $THR $TOTTHR $INPREF $TTYPE $OUTPREF

echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit