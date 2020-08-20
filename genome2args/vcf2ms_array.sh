#!/bin/bash
#$ -N vcf2ms
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=8G

## Specify at submit time, match # of line in GENELIST file
# -t 1-13

echo "_START_$(date)"
# module load Anaconda3/5.3.0

GITPATH='/grid/siepel/home_norepl/mo'

POPFILE=$1 # used for file names
GENELIST=$2 #

GENEL=($(awk '{print $1}' $GENELIST))
CHRL=($(awk '{print $2}' $GENELIST))
FROML=($(awk '{print $3}' $GENELIST))
TOL=($(awk '{print $4}' $GENELIST))

#for ((SGE_TASK_ID=1;SGE_TASK_ID<=13;SGE_TASK_ID++)); do
IDX=$((SGE_TASK_ID-1))
GENE=${GENEL[$IDX]}
CHR=${CHRL[$IDX]}
FROM=${FROML[$IDX]}
TO=${TOL[$IDX]}

echo "Converting ${GENE}_chr${CHR}_${FROM}_${TO} vcf to ms"

${GITPATH}/arg-selection/genome2args/vcf2ms.py vcfs/${POPFILE}_${GENE}_chr${CHR}_${FROM}_${TO}_gt.recode.vcf.gz ms_file/${POPFILE}_${GENE}_chr${CHR}_${FROM}_${TO}_gt $FROM
echo "VCF2HAP_EXITSTAT_$?"

echo "_END_$(date)"
exit