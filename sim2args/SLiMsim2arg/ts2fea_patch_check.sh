#!/bin/bash

## ARRAY TASKS RUNNING ##
PR_JOB_ID=64817910
RUNNING=( $(qstat | grep $PR_JOB_ID | awk '{print $10}') )

## SUCCESSFULLY FINISHED ##
PR_OUTPREF=Trial_fea/trial_swp
RAN=( $(ls ${PR_OUTPREF}_meta_*.npy | awk -F'[_.]' '{print $(NF-1)}') )

EXCLUDE=( "${RUNNING[@]}" "${RAN[@]}" )

echo "EXCLUDE: ${EXCLUDE[@]}"
echo "Total excluded: ${#EXCLUDE[@]}"

# -w, --word-regexp
# Select  only  those  lines containing matches that form whole words.
# -q, --quiet, --silent
# Quiet; do not write anything to standard output.  Exit immediately with zero status if any match is found, even if an error was detected.

for THR in {1..200}
do
    if echo ${EXCLUDE[@]} | grep -q -w $THR; then 
        echo $THR "SKIPPED"
    else
        if [ -f ${PR_OUTPREF}_fea_${THR}.npy ]; then
            ls -lh ${PR_OUTPREF}_fea_${THR}.npy
        fi
        echo $THR "EXECUTED"
    fi
done