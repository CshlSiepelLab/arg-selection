#!/bin/bash

echo "_START_$(date)"

GENELIST="pos_sel_genes.txt"

GENEL=($(awk '{print $1}' $GENELIST))
CHRL=($(awk '{print $2}' $GENELIST))
FROML=($(awk '{print $3}' $GENELIST))
TOL=($(awk '{print $4}' $GENELIST))

for ((SGE_TASK_ID=1;SGE_TASK_ID<=13;SGE_TASK_ID++)); do
	IDX=$((SGE_TASK_ID-1))
	GENE=${GENEL[$IDX]}
	CHR=${CHRL[$IDX]}
	FROM=${FROML[$IDX]}
	TO=${TOL[$IDX]}
	
	#echo Running Relate on ${GENE}_chr${CHR}_${FROM}_${TO}
	./geneTrees_emp.py ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.trees ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}
	echo "_EXITSTAT_$?"

done

echo "_END_$(date)"
