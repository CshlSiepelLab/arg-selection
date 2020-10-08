#!/bin/bash
#$ -N RELATE_array
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=8G

## Specify at submit time, match # of chromosomes
# -t 1-16

echo "_START_$(date)"

GITPATH='/grid/siepel/home_norepl/mo'
RELATE_PATH="/grid/siepel/home_norepl/mo/relate_v1.0.17_x86_64_static"

POP=$1 # e.g. LDR_40
CHR=$SGE_TASK_ID

mkdir ${POP}_chr${CHR}
mv ${POP}_${CHR}.* ${POP}_chr${CHR}/
cp ${POP}.poplabels ${POP}_chr${CHR}/${POP}_${CHR}.poplabels
cp ${POP}.sample ${POP}_chr${CHR}/${POP}_${CHR}.sample

cd ${POP}_chr${CHR}
echo Running Relate on ${POP}_chr${CHR}

# inferred recomb rate 7e-9
${GITPATH}/arg-selection/genome2args/make_uniform_map.py ${POP}_${CHR}.haps 0.7 ${POP}_${CHR}
echo "MAKE_MAP_EXITSTAT_$?"

${RELATE_PATH}/bin/Relate \
	--mode All \
	-m 3e-8 \
	-N 2e8 \
	--haps ${POP}_${CHR}.haps \
	--sample ${POP}_${CHR}.sample \
	--map ${POP}_${CHR}.map \
	-o ${POP}_${CHR}
echo "RELATE_EXITSTAT_$?"

# re-estimate branch lengths
${RELATE_PATH}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
            -i ${POP}_${CHR} \
            -m 3e-8 \
            --poplabels ${POP}_${CHR}.poplabels \
            --threshold 10 \
            -o ${POP}_${CHR}_popsize
echo "POPSIZE_EXITSTAT_$?"

# re-estimate branch length for ENTIRE genealogy
${RELATE_PATH}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
            -i ${POP}_${CHR} \
            -m 3e-8 \
            --poplabels ${POP}_${CHR}.poplabels \
            --threshold 0 \
            --coal ${POP}_${CHR}_popsize.coal \
            --num_iter 1 \
            -o ${POP}_${CHR}_wg
echo "WG_EXITSTAT_$?"

${RELATE_PATH}/bin/RelateFileFormats \
	--mode ConvertToTreeSequence \
	-i ${POP}_${CHR}_wg \
	-o ${POP}_${CHR}_wg
echo "TS_CONVERSION_EXITSTAT_$?"

echo "_END_$(date)"
exit
