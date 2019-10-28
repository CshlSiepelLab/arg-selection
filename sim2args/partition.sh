#!/bin/bash
#$ -N partition_pkl
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=32G


echo "_START_$(date)"

module purge
module load Anaconda3/5.3.0

GITPATH='/sonas-hs/siepel/hpc_norepl/home/mo'
PKL=$1
OUTPREF=$2
THREADS=$3
MODE=$4

# usage: $./partition.py <pkl_path> <out_pref> <no_threads> <mode>
#     - makes directory <out_pref>
#     - outputs <out_pref>_pgv_<thread>.pickle files in directory
#     - calculates and normalizes the iHS score
#     - <mode> can be `s` (sweep) or `n` (neutral)

${GITPATH}/arg-selection/sim2args/partition.py $PKL $OUTPREF $THREADS $MODE

# HANDLEFILE=handles.txt

# HANDLELS=($(awk '{print $1}' $HANDLEFILE))
# THREADLS=($(awk '{print $2}' $HANDLEFILE))

# for LINE in $(seq 0 12); do
# 	echo Partitioning discoal_${HANDLELS[$LINE]}.pkl into ${THREADLS[$LINE]} files.
# 	./partition.py /sonas-hs/siepel/hpc_norepl/home/hijazi/subset/discoal_${HANDLELS[$LINE]}.pkl discoal_subset_${HANDLELS[$LINE]} ${THREADLS[$LINE]}
# 	echo "_EXITSTAT_$?"

# done

echo "_END_$(date)"
exit
