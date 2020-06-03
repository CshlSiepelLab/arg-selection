#!/bin/bash
#$ -N discoal_sim
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_ID_$TASK_ID.o
#$ -e $JOB_ID_$TASK_ID.e
#$ -l m_mem_free=8G

## specify at submit time
# -t 1-3
# -tc 3

# execute `mkdir discoal_${TAG}_${NOCHR}_${COND}` before submitting job

#########################################################
# Simulation of constant Ne model for power analysis    #
#########################################################

echo "_START_$(date)"
DISCOAL_PATH='/sonas-hs/siepel/hpc_norepl/home/mo/discoal-master'

META=$1
NOCHR=$2
TAG=$3
COND=$4 # ("swp" or "neu")

Ne=10000 # constant population size (10000)
pos=0.5 #beneficial mutation positions (0-1)
len=100000 #length of the simulated region
mu=2.5e-8
rho=1.25e-8
SIGMA=40 # time discretization for sweep simulation range: (4, 400),dt=(1/sigma*N)

theta=$(echo $Ne $mu $len | awk '{ printf "%.1f", 4 * $1 * $2 * $3 }')
R=$(echo $Ne $rho $len | awk '{ printf "%.1f", 4 * $1 * $2 * $3 }')

while IFS=$'\t' read -r SIMID SC AF; do
    while :
    do
        if [ $COND == 'swp' ]; then
            sel=$(echo "2*${Ne}*${SC}" | bc -l)
            echo "Attempt: ${DISCOAL_PATH}/discoal $NOCHR 1 $len -t $theta -r $R -N $Ne -i $SIGMA -c $AF -ws 0 -a $sel -x $pos -T > discoal_${TAG}_${NOCHR}_${COND}/discoal_${TAG}_${NOCHR}_${COND}_${SIMID}.discoal"
            ${DISCOAL_PATH}/discoal $NOCHR 1 $len -t $theta -r $R -N $Ne -i $SIGMA -c $AF -ws 0 -a $sel -x $pos -T > discoal_${TAG}_${NOCHR}_${COND}/discoal_${TAG}_${NOCHR}_${COND}_${SIMID}.discoal

        elif [ $COND == 'neu' ]; then
            echo "Attempt: ${DISCOAL_PATH}/discoal $NOCHR 1 $len -t $theta -r $R -T > discoal_${TAG}_${NOCHR}_${COND}/discoal_${TAG}_${NOCHR}_${COND}_${SIMID}.discoal"
            ${DISCOAL_PATH}/discoal $NOCHR 1 $len -t $theta -r $R -T > discoal_${TAG}_${NOCHR}_${COND}/discoal_${TAG}_${NOCHR}_${COND}_${SIMID}.discoal
        fi
        RTCD=$?
        echo "${SIMID}_EXITSTAT_${RTCD}"
        if ((RTCD == 0)); then
            break
        fi
    done
done < "ID_SC_AF/discoal_${META}_${SGE_TASK_ID}_id_sc_af.txt"

echo "_END_$(date)"
exit
