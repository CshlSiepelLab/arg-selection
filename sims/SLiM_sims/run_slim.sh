#!/bin/bash
#$ -N SLiM
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=32G

echo "_START_$(date)"

module purge
module load GCC/7.3.0-2.30

SCRIPT=$2
SEED=$1

echo "RANDOM SEED: $SEED"
echo "RUNNING SCRIPT: $SCRIPT"
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

./slim -s $SEED $SCRIPT

echo "_EXITSTAT_$?"
echo "_END_$(date)"

exit
