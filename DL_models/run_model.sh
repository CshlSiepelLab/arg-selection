#!/bin/bash
#$ -N run_model
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=64G
#$ -l gpu=2

module purge
# module load shared
# module load openmpi-geib-cuda10.2-gcc/3.1.4
# module load tensorflow2-py37-cuda10.2-gcc/2.0.0
module load EBModules
module load OpenMPI/3.1.4-gcccuda-2019b
module load TensorFlow/2.2.0-fosscuda-2019b-Python-3.7.4

GITPATH='/grid/siepel/home_norepl/mo'
echo "_START_$(date)"

HAPS=$1
FEA_PREF=$2
NO_THR=$3
MODEL=$4
OUT_PATH=$5

${GITPATH}/arg-selection/DL_models/SIA_class_pred.py $HAPS $FEA_PREF $NO_THR $MODEL $OUT_PATH

echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit