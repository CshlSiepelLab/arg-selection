#!/bin/bash
#$ -N SLiM_array
#$ -S /bin/bash
#$ -cwd
#$ -t 1-20
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=32G

echo "_START_$(date)"

module purge
module load GCC/7.3.0-2.30

#mapfile -t SAMPIDS < $1
#INDV=${SAMPIDS[$((SGE_TASK_ID-1))]}
INDV=$1
SUFFIX=$2
MUTGEN=$3
SCMIN=$4
SCMAX=$5
MINAF=$6
SCRIPT=$7
RUNS=$8
LASTIDX=$9

echo "PARAMETER SET: $INDV"
echo "SLIM SCRIPT: $SCRIPT"
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
	./slim -s $RANDOM -t -d "paramF='slim_params/${INDV}.psmc_scaled.param'" -d "outPref='slim_output/${INDV}_${SUFFIX}/hdswp_${INDV}.psmc_$((LASTIDX+(SGE_TASK_ID-1)*RUNS+sim))'" -d "mutgen=$MUTGEN" -d "sc_min=$SCMIN" -d "sc_max=$SCMAX" -d "min_AF=$MINAF" $SCRIPT

done

echo "_EXITSTAT_$?"
echo "_END_$(date)"

exit
