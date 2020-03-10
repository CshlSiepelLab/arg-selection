#!/bin/bash
#$ -N SLiM_patch
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=20G

echo "_START_$(date)"

module purge
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load GSL/2.5

SLIMDIR="/sonas-hs/siepel/hpc_norepl/home/mo/maiz/build"
GITPATH="/sonas-hs/siepel/hpc_norepl/home/mo/arg-selection"
SCRIPT="${GITPATH}/sims/SLiM_sims/sweep_treeseq.slim"

PARAMF=$1
NOCHR=$2
minSC=$3
maxSC=$4
minAF=$5
maxAF=$6
OUTPREF=$7
START=$8
END=$9

#usage: slim -v[ersion] | -u[sage] | -testEidos | -testSLiM |
#   [-l[ong]] [-s[eed] <seed>] [-t[ime]] [-m[em]] [-M[emhist]] [-x]
#   [-d[efine] <def>] [<script file>]

#   -l[ong]          : long (i.e.) verbose output (format may change)
#   -s[eed] <seed>   : supply an initial random number seed for SLiM
#   -t[ime]          : print SLiM's total execution time (in user clock time)
#   -m[em]           : print SLiM's peak memory usage
#   -M[emhist]       : print a histogram of SLiM's memory usage
#   -x               : disable SLiM's runtime safety/consistency checks
#   -d[efine] <def>  : define an Eidos constant, such as "mu=1e-7"
#   <script file>    : the input script file (stdin may be used instead)

echo Patching $START to $END

for RUN_ID in $(seq $START $END); do
	while :
	do
		${SLIMDIR}/slim -s $(tr -cd "[:digit:]" < /dev/urandom | head -c 10) -d "paramF='${PARAMF}'" -d "outPref='${OUTPREF}_${RUN_ID}'" -d "sc_min=${minSC}" -d "sc_max=${maxSC}" -d "min_AF=${minAF}" -d "max_AF=${maxAF}" $SCRIPT
		SLIM_RTCD=$?
		echo "${RUN_ID}_SLiM_EXITSTAT_${SLIM_RTCD}"
		if ((SLIM_RTCD != 0)); then continue ; fi
		#/usr/bin/time -f "RSS=%M elapsed=%E"
		${GITPATH}/sims/SLiM_sims/recapitation.py ${OUTPREF}_${RUN_ID}.trees ${PARAMF} $NOCHR ${OUTPREF}_${RUN_ID}_samp
		REC_RTCD=$?
		echo "${RUN_ID}_recap_EXITSTAT_${REC_RTCD}"
		if ((REC_RTCD == 0)); then break ; fi
	done
	# delete the intermediate tree file (takes up too much space)
	rm ${OUTPREF}_${RUN_ID}.trees
done

echo "_END_$(date)"
exit
