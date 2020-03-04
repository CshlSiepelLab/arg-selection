#!/bin/bash
#$ -N SLiM_array
#$ -S /bin/bash
#$ -cwd
#$ -t 1-100
#$ -tc 50
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
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
LASTIDX=$8
RUNS=$9 # no of new runs PER THREAD

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

for sim in $(seq 1 $RUNS); do
	RUN_ID=$((LASTIDX+(SGE_TASK_ID-1)*RUNS+sim))
	${SLIMDIR}/slim -s $(tr -cd "[:digit:]" < /dev/urandom | head -c 10) -t -m -d "paramF='${PARAMF}'" -d "outPref='${OUTPREF}_${RUN_ID}'" -d "sc_min=${minSC}" -d "sc_max=${maxSC}" -d "min_AF=${minAF}" -d "max_AF=${maxAF}" $SCRIPT
	echo "${RUN_ID}_SLiM_EXITSTAT_$?"
	/usr/bin/time -f "RSS=%M elapsed=%E" ${GITPATH}/sims/SLiM_sims/recapitation.py ${OUTPREF}_${RUN_ID}.trees ${PARAMF} $NOCHR ${OUTPREF}_${RUN_ID}_samp
	echo "${RUN_ID}_recap_EXITSTAT_$?"
	# delete the intermediate tree file (takes up too much space)
	rm ${OUTPREF}_${RUN_ID}.trees
done

echo "_END_$(date)"
exit
