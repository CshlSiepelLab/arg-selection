#!/bin/bash
#$ -N hmp2haps_samp
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=8G

echo "_START_$(date)"
#module load Anaconda3/5.3.0

GITPATH='/grid/siepel/home_norepl/mo'

HAPMAP_PATH=$1
OUTPREF=$2

${GITPATH}/arg-selection/genome2args/hmp2hap.py $HAPMAP_PATH $OUTPREF
echo "HMP2HAP_EXITSTAT_$?"

${GITPATH}/arg-selection/genome2args/make_poplabels.py ${OUTPREF}.sample CEU EUR $OUTPREF
echo "MAKE_POPLABELS_EXITSTAT_$?"

echo "_END_$(date)"
exit