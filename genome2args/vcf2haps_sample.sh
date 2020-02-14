GITPATH='/sonas-hs/siepel/hpc_norepl/home/mo'
RELATE_PATH="/sonas-hs/siepel/hpc_norepl/home/mo/relate_v1.0.17_x86_64_static"

module load Anaconda3/5.3.0

while IFS=$'\t' read -r GENE CHR FROM TO; do
	mkdir $GENE
    echo Converting: vcfs/CEU_99.txt_${GENE}_chr${CHR}_${FROM}_${TO}_gt.recode.vcf.gz > ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log
    ${RELATE_PATH}/bin/RelateFileFormats \
	--mode ConvertFromVcf \
	--haps ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}_dummy.haps \
	--sample ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.sample \
	-i vcfs/CEU_99.txt_${GENE}_chr${CHR}_${FROM}_${TO}_gt.recode \
	--chr $CHR >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log
    echo "RELATEFF_EXITSTAT_$?" >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log

	${GITPATH}/arg-selection/genome2args/vcf2hap.py vcfs/CEU_99.txt_${GENE}_chr${CHR}_${FROM}_${TO}_gt.recode.vcf.gz ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.haps >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log
	echo "VCF2HAP_EXITSTAT_$?" >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log

	${GITPATH}/arg-selection/genome2args/make_poplabels.py ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.sample CEU EUR ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}
	echo "MAKE_POPLABELS_EXITSTAT_$?" >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log

	${GITPATH}/arg-selection/genome2args/make_uniform_map.py ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.haps 1.25 ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}
	echo "MAKE_MAP_EXITSTAT_$?" >> ${GENE}/${GENE}_chr${CHR}_${FROM}_${TO}.log

done < ${GITPATH}/arg-selection/genome2args/pos_sel_genes.txt