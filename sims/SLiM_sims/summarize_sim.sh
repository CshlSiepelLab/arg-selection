OUTPREF=$1
JOBID=$2
TOTTHR=$3

echo -e 'THR\tSLIM_SUC\tSLIM_FAIL\tRECAP_SUC\tRECAP_FAIL\tMETA_LINES'
for ((t=1;t<=TOTTHR;t++)); do
    SLIM_SUC=$(grep -c SLiM_EXITSTAT_0 ${JOBID}_${t}.o)
    SLIM_FAL=$(grep -c SLiM_EXITSTAT_1 ${JOBID}_${t}.o)
    RECAP_SUC=$(grep -c recap_EXITSTAT_0 ${JOBID}_${t}.o)
    RECAP_FAIL=$(grep -c recap_EXITSTAT_1 ${JOBID}_${t}.o)
    META_LINS=$(grep -c %% ${JOBID}_${t}.o)
    echo -e $t'\t'$SLIM_SUC'\t'$SLIM_FAL'\t'$RECAP_SUC'\t'$RECAP_FAIL'\t'$META_LINS
done
