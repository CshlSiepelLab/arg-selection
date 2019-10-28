#!/usr/bin/env python3

import sys, os
import pickle
import numpy as np
#import pandas as pd
import allel

helpMsg = '''
        usage: $./partition.py <pkl_path> <out_pref> <no_threads> <mode>
            - makes directory <out_pref>
            - outputs <out_pref>_pgv_<thread>.pickle files in directory
            - calculates and normalizes the iHS score
            - <mode> can be `s` (sweep) or `n` (neutral)
'''

def all_iHS(gt_mtx, pos):

    ppos = np.around(pos * 100000) # convert position to coordinate in 100kb region
    hArr = allel.HaplotypeArray(gt_mtx)
    acArr = hArr.count_alleles() 

    biall05 = acArr.is_biallelic_01(min_mac=np.ceil(0.05*gt_mtx.shape[1]))
    IHS = allel.ihs(hArr, ppos, min_maf=0.05, include_edges=True) # TBD

    varStats = np.stack((pos[biall05], IHS[biall05], acArr[:, 1][biall05]), axis=-1)

    return varStats

def main(args):
    if len(args) != 5:    #4 arguments
        return helpMsg

    pkl_path = args[1]
    out_pref = args[2]
    no_threads = int(args[3])
    mode = args[4]

    if mode == 's':    
        with open(pkl_path, 'rb') as f:
            list_sel_coef, list_freq, list_variant, list_tree, list_geno, list_pos = pickle.load(f)
    elif mode == 'n':
        with open(pkl_path, 'rb') as f:
            list_sel_coef, list_freq, list_variant, list_tree, list_geno, list_pos, list_variant_pos = pickle.load(f)

    ## iHS calculation ##
    no_sims = len(list_pos) # number of simulations
    sim_ids = np.arange(no_sims)

    iHS_df = np.empty((0, 4))
    for i in range(no_sims):
        r_idx = sim_ids[i]
        iHS_var = all_iHS(list_geno[r_idx], list_pos[r_idx])
        iHS_df = np.vstack((iHS_df, np.hstack((np.repeat(r_idx, iHS_var.shape[0]).reshape(iHS_var.shape[0],1), iHS_var))))

    std_iHS, _ = allel.standardize_by_allele_count(iHS_df[:, 2], 198-iHS_df[:, 3], n_bins=50)
    iHS_df[:, 2] = std_iHS

    #df_cols=["ID", "pos", "iHS", "aac"]
    #iHS_pd_df = pd.DataFrame(iHS_df, columns=df_cols)

    # Partition
    tasks = len(list_pos)//no_threads
    os.mkdir(out_pref, 0o755)

    for thread in range(1, no_threads+1):
        a_idx = (thread-1)*tasks
        if thread == no_threads:
            b_idx = len(list_pos)
        else:
            b_idx = thread*tasks

        iHS_sub = iHS_df[iHS_df[:, 0]>=a_idx & iHS_df[:, 0]<b_idx, :]

        if mode == 's':
            with open(out_pref+'/'+out_pref+'_pgv_'+str(thread)+'.pickle', 'wb') as f:
                pickle.dump((list_pos[a_idx:b_idx], list_geno[a_idx:b_idx], list_variant[a_idx:b_idx], iHS_sub), f, pickle.HIGHEST_PROTOCOL)
        elif mode == 'n':
            with open(out_pref+'/'+out_pref+'_pgv_'+str(thread)+'.pickle', 'wb') as f:
                pickle.dump((list_pos[a_idx:b_idx], list_geno[a_idx:b_idx], list_variant[a_idx:b_idx], list_variant_pos[a_idx:b_idx], iHS_sub), f, pickle.HIGHEST_PROTOCOL)
    # TO BE TESTED
    return 0

sys.exit(main(sys.argv))