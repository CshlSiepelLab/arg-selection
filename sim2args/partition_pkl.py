#!/usr/bin/env python3

import sys, os
import pickle
import numpy as np
#import pandas as pd
#import allel

helpMsg = '''
        usage: $./partition_pkl.py <pkl_path> <pkl_pref> <part_i> <no_out_p_in>
            - need to manually make directory <pkl_pref> before submitting the job
            - <part_i>: thread # of input partition
            - <no_out_p_in>: # of output partitions per input
            - outputs <pkl_pref>_<thread>.pkl files in respective directory
'''

def main(args):
    if len(args) != 5:    #4 arguments
        return helpMsg

    path = args[1]
    pkl_pref = args[2]
    partI = int(args[3])
    no_partO = int(args[4])

    with open(path+'/'+pkl_pref+'_'+str(partI)+'.pkl', 'rb') as f:
    	list_sel_coef, list_freq, list_variant, _, list_geno, list_pos, list_variant_pos, _, theta, rho, Ne, _ = pickle.load(f)

    # Partition
    no_sims = len(list_sel_coef)
    tasks = no_sims//no_partO
    #os.mkdir(pkl_pref, 0o755) # <= make this directory beforehand

    print("%% Splitting:", pkl_pref+'_'+str(partI)+'.pkl', no_sims, "to", no_partO, flush=True)
    for thread in range(no_partO):
        a_idx = thread*tasks
        if thread == no_partO-1:
            b_idx = no_sims
        else:
            b_idx = (thread+1)*tasks

        glob_thr = partI*no_partO+thread
        print("%%", a_idx, b_idx, "to global thread", glob_thr, flush=True)
        with open(pkl_pref+'/'+pkl_pref+'_'+str(glob_thr)+'.pkl', 'wb') as f:
            pickle.dump((list_sel_coef[a_idx:b_idx], list_freq[a_idx:b_idx], list_pos[a_idx:b_idx], list_geno[a_idx:b_idx], list_variant[a_idx:b_idx], list_variant_pos[a_idx:b_idx], theta[a_idx:b_idx], rho[a_idx:b_idx], Ne[a_idx:b_idx]), f, pickle.HIGHEST_PROTOCOL)

    return 0

sys.exit(main(sys.argv))
