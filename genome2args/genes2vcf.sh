POP=$1 # file name
GENELIST=$2 # file name

GITPATH='/sonas-hs/siepel/hpc_norepl/home/mo'
OUTPATH='/sonas-hs/siepel/hpc_norepl/home/mo/SIA_GBR_empirical/vcfs' # change this manually

# vcf format uses 1-based indexing
# IFS - Internal Field Separator

while IFS=$'\t' read -r GENE CHR FROM TO; do
	echo executing: nohup ${GITPATH}/arg-selection/genome2args/xtract_local_vcf.sh $POP $CHR $FROM $TO $GENE $OUTPATH '&'
	nohup ${GITPATH}/arg-selection/genome2args/xtract_local_vcf.sh $POP $CHR $FROM $TO $GENE $OUTPATH &
	sleep 1 
done < "$GENELIST"