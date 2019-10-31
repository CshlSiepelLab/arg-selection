module load Anaconda3/5.3.0

while IFS=$'\t' read -r GENE CHR FROM TO; do
        echo Converting: vcfs/CEU_99.txt_${GENE}_chr${CHR}_${FROM}_${TO}_gt.recode.vcf.gz > ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log
	#mkdir $GENE
        relate_v1.0.16_x86_64_static/bin/RelateFileFormats \
		--mode ConvertFromVcf \
		--haps ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}_dummy.haps \
		--sample ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.sample \
		-i vcfs/CEU_99.txt_${GENE}_chr${CHR}_${FROM}_${TO}_gt.recode \
		--chr $CHR >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log
        echo "_EXITSTAT_$?" >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log
	
	./vcf2hap.py vcfs/CEU_99.txt_${GENE}_chr${CHR}_${FROM}_${TO}_gt.recode.vcf.gz ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.haps >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log
	echo "_EXITSTAT_$?" >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log

done < pos_sel_genes.txt

