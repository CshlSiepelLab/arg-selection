#!/bin/bash
#$ -N slimSims2ts
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=32G

echo "_START_$(date)"

module purge
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load GSL/2.5

GITPATH="/sonas-hs/siepel/hpc_norepl/home/mo/arg-selection"

PARAMF=$1
SCALE=$2
MODE=$3 # <`n`/`s`>
META=$4 # <no_sims/meta_file_path>
THR=$5
TOTTHR=$6
INPREF=$7
OUTPREF=$8

#echo Patching $START to $END

${GITPATH}/sim2args/SLiMsim2arg/slimSims2ts.py $PARAMF $SCALE $MODE $META $THR $TOTTHR $INPREF $OUTPREF

echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit