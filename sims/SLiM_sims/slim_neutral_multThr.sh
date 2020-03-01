#!/bin/bash
#$ -N SLiM_array
#$ -S /bin/bash
#$ -cwd
#$ -t 1-5
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=32G

echo "_START_$(date)"

module purge
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load GSL/2.5

SLIMDIR="/sonas-hs/siepel/hpc_norepl/home/mo/maiz/build"
GITPATH="/sonas-hs/siepel/hpc_norepl/home/mo/arg-selection"
SCRIPT="${GITPATH}/sims/SLiM_sims/neutral.slim"

PARAMF=$1
NOCHR=$2
OUTPREF=$3
LASTIDX=$4
RUNS=$5 # no of new runs PER THREAD

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
	${SLIMDIR}/slim -s $RANDOM -t -m -d "paramF='${PARAMF}'" -d "outPref='${OUTPREF}_${RUN_ID}_temp'" $SCRIPT
	echo "${RUN_ID}_SLiM_EXITSTAT_$?"
	${GITPATH}/sims/SLiM_sims/recapitation.py ${OUTPREF}_${RUN_ID}_temp.trees ${PARAMF} $NOCHR ${OUTPREF}_${RUN_ID}_samp
	echo "${RUN_ID}_recap_EXITSTAT_$?"
	# delete the tree file
	rm ${OUTPREF}_${RUN_ID}_temp.trees
done

echo "_END_$(date)"
exit