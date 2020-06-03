#!/usr/bin/env python3

import sys, os
import pickle
import numpy as np
#import pandas as pd
#import allel

helpMsg = '''
        usage: $./partition.py <pkl_path> <pkl_pref> <part_i> <no_out_p_in>
            - need to manually make directory <pkl_pref> before submitting the job
            - <part_i>: # of input partitions
            - <no_out_p_in>: # of output partitions per input
            - outputs <out_pref>_pgv_<thread>.pkl files in respective directory
'''

def main(args):
    if len(args) != 5:    #4 arguments
        return helpMsg

    path = args[1]
    pkl_pref = args[2]
    partI = int(args[3])
    no_partO = int(args[4])

    #for p_in in range(partI):
        # pickle fmt: list_sel_coef(0, ls), list_freq(1, ls), list_variant(2, nparr), 
        # trees_of_interest(3, ls), list_geno(4, ls), list_pos(5, ls), list_variant_pos(6, ls), length_ARG(7, ls)
    with open(path+'_'+str(partI)+'.pkl', 'rb') as f:
    	_, _, list_variant, _, list_geno, list_pos, list_variant_pos, _, list_mu, list_rho, _ = pickle.load(f) # sample vary in mu & rho

    # Partition
    no_sims = len(list_pos)
    tasks = no_sims//no_partO
    #os.mkdir(pkl_pref, 0o755) # <= make this directory beforehand

    print("%% Splitting:", pkl_pref+'_'+str(partI)+'.pkl', no_sims, "to", no_partO, flush=True)
    for thread in range(no_partO):
        a_idx = thread*tasks
        if thread == no_partO-1:
            b_idx = no_sims
        else:
            b_idx = (thread+1)*tasks

        #iHS_sub = siHS_df[(siHS_df[:, 0]>=a_idx) & (siHS_df[:, 0]<b_idx), :]
        iHS_sub = None
        glob_thr = partI*no_partO+thread
        print("%%", a_idx, b_idx, "to global thread", glob_thr, flush=True)
        with open(pkl_pref+'/'+pkl_pref+'_pgv_'+str(glob_thr)+'.pkl', 'wb') as f:
            #pickle.dump((list_pos[a_idx:b_idx], list_geno[a_idx:b_idx], list_variant[a_idx:b_idx], list_variant_pos[a_idx:b_idx], list_mu[a_idx:b_idx], list_rho[a_idx:b_idx], list_Ne[a_idx:b_idx], iHS_sub), f, pickle.HIGHEST_PROTOCOL)
            pickle.dump((list_pos[a_idx:b_idx], list_geno[a_idx:b_idx], list_variant[a_idx:b_idx], list_variant_pos[a_idx:b_idx], list_mu[a_idx:b_idx], list_rho[a_idx:b_idx], iHS_sub), f, pickle.HIGHEST_PROTOCOL)

    return 0

sys.exit(main(sys.argv))
