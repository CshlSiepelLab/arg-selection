#!/bin/bash
#$ -N RELATE_array
#$ -S /bin/bash
#$ -cwd
#$ -t 1-13
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=16G

echo "_START_$(date)"

GENELIST="pos_sel_genes.txt"

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

cd $GENE
echo Running Relate on ${GENE}_chr${CHR}_${FROM}_${TO}

#../relate_v1.0.16_x86_64_static/scripts/PrepareInputFiles/PrepareInputFiles.sh \
#	--haps ${GENE}_chr${CHR}_${FROM}_${TO}.haps \
#	--sample ${GENE}_chr${CHR}_${FROM}_${TO}.sample \
#	--ancestor ../ancestral_fasta/homo_sapiens_ancestor_${CHR}.fa \
#	-o ${GENE}_chr${CHR}_${FROM}_${TO}_input

#gunzip ${GENE}_chr${CHR}_${FROM}_${TO}_input.haps.gz
#gunzip ${GENE}_chr${CHR}_${FROM}_${TO}_input.sample.gz

../make_uniform_map.py ${GENE}_chr${CHR}_${FROM}_${TO}.haps 1.25 ${GENE}_chr${CHR}_${FROM}_${TO}
echo "MAKE_MAP_EXITSTAT_$?"

../relate_v1.0.16_x86_64_static/bin/Relate \
	--mode All \
	-m 2.5e-8 \
	-N 376176 \
	--haps ${GENE}_chr${CHR}_${FROM}_${TO}.haps \
	--sample ${GENE}_chr${CHR}_${FROM}_${TO}.sample \
	--map ${GENE}_chr${CHR}_${FROM}_${TO}.map \
	-o ${GENE}_chr${CHR}_${FROM}_${TO}
echo "RELATE_EXITSTAT_$?"

../relate_v1.0.16_x86_64_static/bin/RelateFileFormats \
	--mode ConvertToTreeSequence \
	-i ${GENE}_chr${CHR}_${FROM}_${TO} \
	-o ${GENE}_chr${CHR}_${FROM}_${TO}

echo "TS_CONVERSION_EXITSTAT_$?"
echo "_END_$(date)"

exit
