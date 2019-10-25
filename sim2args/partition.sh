#!/bin/bash
#$ -N partition_pkl
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=32G


echo "_START_$(date)"

module purge
module load Anaconda3/5.3.0

#HANDLE=$1
#THREADS=$2

#        usage: $./partition.py <pkl_path> <out_pref> <no_threads>
#            - makes directory <out_pref>
#            - outputs <out_pref>_pgv_<thread>.pickle files in directory

HANDLEFILE=handles.txt

HANDLELS=($(awk '{print $1}' $HANDLEFILE))
THREADLS=($(awk '{print $2}' $HANDLEFILE))

for LINE in $(seq 0 12); do
	echo Partitioning discoal_${HANDLELS[$LINE]}.pkl into ${THREADLS[$LINE]} files.
	./partition.py /sonas-hs/siepel/hpc_norepl/home/hijazi/subset/discoal_${HANDLELS[$LINE]}.pkl discoal_subset_${HANDLELS[$LINE]} ${THREADLS[$LINE]}
	echo "_EXITSTAT_$?"

done

echo "_END_$(date)"
exit
