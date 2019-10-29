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
FILEPATH=$1 #path to pickle files, ending with '/''
SWPPREF=$2
NEUPREF=$3
THREADS=$4

# usage: $./partition.py <pkl_path> <swp_pkl_pref> <neu_pkl_pref> <no_threads>
#     - makes directories <swp_pkl_pref> <neu_pkl_pref>
#     - outputs <out_pref>_pgv_<thread>.pkl files in respective directory
#     - calculates and normalizes the iHS score
#     - sweep and neutral files needed simultaneously for normalization purpose

${GITPATH}/arg-selection/sim2args/partition.py $FILEPATH $SWPPREF $NEUPREF $THREADS
echo "_EXITSTAT_$?"

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
