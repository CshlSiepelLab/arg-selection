#!/usr/bin/env python3

import sys, os
import pickle
import numpy as np
#import pandas as pd
import allel

helpMsg = '''
        usage: $./partition.py <pkl_path> <swp_pkl_pref> <neu_pkl_pref> <part_i> <part_o>
            - makes directories <swp_pkl_pref> <neu_pkl_pref>
            - <part_i>: # of input partitions
            - <part_o>: # of output partitions
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
    if len(args) != 6:    #5 arguments
        return helpMsg

    path = args[1]
    swp_pref = args[2]
    neu_pref = args[3]
    partI = int(args[4])
    partO = int(args[5])
    #mode = args[4]

    slist_pos = []
    slist_geno = []
    slist_variant = np.empty((0, 198)) # of haplotypes = 198
    slist_variant_pos = []
    nlist_pos = []
    nlist_geno = []
    nlist_variant = np.empty((0, 198)) # of haplotypes = 198
    nlist_variant_pos = []

    for p_in in range(partI):
        # pickle fmt: list_sel_coef(0, ls), list_freq(1, ls), list_variant(2, nparr), 
        # trees_of_interest(3, ls), list_geno(4, ls), list_pos(5, ls), list_variant_pos(6, ls), length_ARG(7, ls)
        with open(path+swp_pref+'_'+str(p_in)+'.pkl', 'rb') as f:
            sllist = pickle.load(f)
            slist_pos += sllist[5]
            slist_geno += sllist[4]
            slist_variant = np.vstack((slist_variant, sllist[2]))
            slist_variant_pos += sllist[6]

        with open(path+neu_pref+'_'+str(p_in)+'.pkl', 'rb') as f:
            nllist = pickle.load(f)
            nlist_pos += nllist[5]
            nlist_geno += nllist[4]
            nlist_variant = np.vstack((nlist_variant, nllist[2]))
            nlist_variant_pos += nllist[6]

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
    tasks = no_ssims//partO
    os.mkdir(swp_pref, 0o755)

    for thread in range(1, partO+1):
        a_idx = (thread-1)*tasks
        if thread == partO:
            b_idx = no_ssims
        else:
            b_idx = thread*tasks

        iHS_sub = siHS_df[(siHS_df[:, 0]>=a_idx) & (siHS_df[:, 0]<b_idx), :]
        with open(swp_pref+'/'+swp_pref+'_pgv_'+str(thread)+'.pkl', 'wb') as f:
            pickle.dump((slist_pos[a_idx:b_idx], slist_geno[a_idx:b_idx], slist_variant[a_idx:b_idx], slist_variant_pos[a_idx:b_idx], iHS_sub), f, pickle.HIGHEST_PROTOCOL)

    tasks = no_nsims//partO
    os.mkdir(neu_pref, 0o755)

    for thread in range(1, partO+1):
        a_idx = (thread-1)*tasks
        if thread == partO:
            b_idx = no_nsims
        else:
            b_idx = thread*tasks

        iHS_sub = niHS_df[(niHS_df[:, 0]>=a_idx) & (niHS_df[:, 0]<b_idx), :]
        with open(neu_pref+'/'+neu_pref+'_pgv_'+str(thread)+'.pkl', 'wb') as f:
            pickle.dump((nlist_pos[a_idx:b_idx], nlist_geno[a_idx:b_idx], nlist_variant[a_idx:b_idx], nlist_variant_pos[a_idx:b_idx], iHS_sub), f, pickle.HIGHEST_PROTOCOL)

    return 0

sys.exit(main(sys.argv))
