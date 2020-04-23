#!/bin/bash
#$ -N SLiM_array
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=64G

# -t 1-50
# -tc 50

echo "_START_$(date)"

module purge
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load GSL/2.5

SLIMDIR="/sonas-hs/siepel/hpc_norepl/home/mo/maiz/build"
GITPATH="/sonas-hs/siepel/hpc_norepl/home/mo/arg-selection"
SCRIPT="${GITPATH}/sims/SLiM_sims/sweep_test.slim"

PARAMF=$1
NOCHR=$2
SCF=$3 # file containing sel coef to be sampled, one per line
minAF=$4
maxAF=$5
HNDL=$6
OUTPREF=${HNDL}/${HNDL}
LASTIDX=$7
RUNS=$8 # no of new runs PER THREAD

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
    if [ -f ${OUTPREF}_${RUN_ID}_samp.trees ]; then
        echo "$RUN_ID was successful, SKIPPING"
        continue
    fi
    SELCOEF=$(shuf -n 1 $SCF)
    while :
    do
        ${SLIMDIR}/slim -s $(tr -cd "[:digit:]" < /dev/urandom | head -c 10) -d "paramF='${PARAMF}'" -d "outPref='${OUTPREF}_${RUN_ID}_temp'" -d "SC=${SELCOEF}" -d "min_AF=${minAF}" -d "max_AF=${maxAF}" $SCRIPT | tee ${JOB_ID}_${SGE_TASK_ID}.buf
        SLIM_RTCD=${PIPESTATUS[0]}
        echo "${RUN_ID}_SLiM_EXITSTAT_${SLIM_RTCD}"
        if ((SLIM_RTCD != 0)); then continue ; fi
        #if ((SLIM_RTCD == 0)); then break ; fi
        ATTEMPTS=0
        while :
        do
            #/usr/bin/time -f "RSS=%M elapsed=%E"
            ${GITPATH}/sims/SLiM_sims/recapitation.py ${OUTPREF}_${RUN_ID}_temp.trees ${PARAMF} $NOCHR ${OUTPREF}_${RUN_ID}_samp
            REC_RTCD=$?
            echo "${RUN_ID}_recap_EXITSTAT_${REC_RTCD}"
            if ((REC_RTCD == 0)); then
                grep "%%" ${JOB_ID}_${SGE_TASK_ID}.buf >> ${HNDL}_${SGE_TASK_ID}.meta
                echo "${RUN_ID}_SUCCESS" >> ${HNDL}_${SGE_TASK_ID}.0exit
                break 2
            fi
            echo "Recap attempt:${ATTEMPTS} FAILED"
            ((ATTEMPTS=ATTEMPTS+1))
            if ((ATTEMPTS > 20)); then break ; fi
        done
    done
    # delete the intermediate tree file (takes up too much space)
    rm ${OUTPREF}_${RUN_ID}_temp.trees
done

echo "_END_$(date)"
exit
