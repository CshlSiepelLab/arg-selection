#!/bin/bash
#$ -N SLiM_array
#$ -S /bin/bash
#$ -cwd
#$ -t 1-100
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=4G

echo "_START_$(date)"

module purge
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load GSL/2.5

SLIMDIR="/sonas-hs/siepel/hpc_norepl/home/mo/maiz/build"
GITPATH="/sonas-hs/siepel/hpc_norepl/home/mo/arg-selection"
SCRIPT="${GITPATH}/sims/SLiM_sims/sweep_survey.slim"

PARAMF=$1
RUNS=$2 # no of new runs PER THREAD
MUTGENEL=$3 # ealiest gen to introduce mutation
MUTGENLT=$4 # latest gen to introduce mutation
SCMIN=$5 # min SCALED sel. coef
SCMAX=$6 # max SCALED sel. coef

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
	${SLIMDIR}/slim -s $RANDOM -t -m -d "paramF='${PARAMF}'" -d "mutgen_early=$MUTGENEL" -d "mutgen_late=$MUTGENLT" -d "sc_min=$SCMIN" -d "sc_max=$SCMAX" $SCRIPT
	echo "_EXITSTAT_$?"
done

echo "_END_$(date)"
exit
