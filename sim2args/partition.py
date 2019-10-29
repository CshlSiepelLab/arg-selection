#!/usr/bin/env python3

import sys, os
import pickle
import numpy as np
#import pandas as pd
import allel

helpMsg = '''
        usage: $./partition.py <pkl_path> <swp_pkl_pref> <neu_pkl_pref> <no_threads>
            - makes directories <swp_pkl_pref> <neu_pkl_pref>
            - outputs <out_pref>_pgv_<thread>.pkl files in respective directory
            - calculates and normalizes the iHS score
            - sweep and neutral files needed simultaneously for normalization purpose
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

    path = args[1]
    swp_pref = args[2]
    neu_pref = args[3]
    no_threads = int(args[4])
    #mode = args[4]

    with open(path+swp_pref+'.pkl', 'rb') as f:
        slist_sel_coef, slist_freq, slist_variant, slist_tree, slist_geno, slist_pos, slist_variant_pos = pickle.load(f)

    with open(path+neu_pref+'.pkl', 'rb') as f:
        nlist_sel_coef, nlist_freq, nlist_variant, nlist_tree, nlist_geno, nlist_pos, nlist_variant_pos = pickle.load(f)

    ## iHS calculation ##
    no_ssims = len(slist_pos) # number of simulations
    ssim_ids = np.arange(no_ssims)

    siHS_df = np.empty((0, 4))
    for i in range(no_ssims):
        r_idx = ssim_ids[i]
        iHS_var = all_iHS(slist_geno[r_idx], slist_pos[r_idx])
        siHS_df = np.vstack((siHS_df, np.hstack((np.repeat(r_idx, iHS_var.shape[0]).reshape(iHS_var.shape[0],1), iHS_var))))

    no_nsims = len(nlist_pos) # number of simulations
    nsim_ids = np.arange(no_nsims)

    niHS_df = np.empty((0, 4))
    for i in range(no_nsims):
        r_idx = nsim_ids[i]
        iHS_var = all_iHS(nlist_geno[r_idx], nlist_pos[r_idx])
        niHS_df = np.vstack((niHS_df, np.hstack((np.repeat(r_idx, iHS_var.shape[0]).reshape(iHS_var.shape[0],1), iHS_var))))

    std_iHS, _ = allel.standardize_by_allele_count(np.concatenate((siHS_df[:, 2], niHS_df[:, 2])), 198-np.concatenate((siHS_df[:, 3], niHS_df[:, 3])), n_bins=50, diagnostics=False)
    siHS_df[:, 2] = std_iHS[:siHS_df.shape[0]]
    niHS_df[:, 2] = std_iHS[siHS_df.shape[0]:]

    #df_cols=["ID", "pos", "iHS", "aac"]
    #iHS_pd_df = pd.DataFrame(iHS_df, columns=df_cols)

    # Partition
    tasks = no_ssims//no_threads
    os.mkdir(swp_pref, 0o755)

    for thread in range(1, no_threads+1):
        a_idx = (thread-1)*tasks
        if thread == no_threads:
            b_idx = no_ssims
        else:
            b_idx = thread*tasks

        iHS_sub = siHS_df[(siHS_df[:, 0]>=a_idx) & (siHS_df[:, 0]<b_idx), :]
        with open(swp_pref+'/'+swp_pref+'_pgv_'+str(thread)+'.pkl', 'wb') as f:
            pickle.dump((slist_pos[a_idx:b_idx], slist_geno[a_idx:b_idx], slist_variant[a_idx:b_idx], slist_variant_pos[a_idx:b_idx], iHS_sub), f, pickle.HIGHEST_PROTOCOL)

    tasks = no_nsims//no_threads
    os.mkdir(neu_pref, 0o755)

    for thread in range(1, no_threads+1):
        a_idx = (thread-1)*tasks
        if thread == no_threads:
            b_idx = no_nsims
        else:
            b_idx = thread*tasks

        iHS_sub = niHS_df[(niHS_df[:, 0]>=a_idx) & (niHS_df[:, 0]<b_idx), :]
        with open(neu_pref+'/'+neu_pref+'_pgv_'+str(thread)+'.pkl', 'wb') as f:
            pickle.dump((nlist_pos[a_idx:b_idx], nlist_geno[a_idx:b_idx], nlist_variant[a_idx:b_idx], nlist_variant_pos[a_idx:b_idx], iHS_sub), f, pickle.HIGHEST_PROTOCOL)

    return 0

sys.exit(main(sys.argv))
