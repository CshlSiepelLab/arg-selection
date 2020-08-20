#!/usr/bin/env python3

import sys, os
import pickle
import numpy as np
#import pandas as pd
#import allel

helpMsg = '''
        usage: $./merge_pkl.py <pkl_path> <pkl_pref> <no_part_i> <no_in_p_out>
            - need to manually make directory <pkl_pref> before submitting the job
            - <part_i>: tot # of input partition to merge per output
            - outputs <pkl_pref>_<thread>.pkl files in respective directory
'''

# bird - 24; CEU - 198; GBR - 182
CHRS=24 # number of chromosomes in simulation


def main(args):
    if len(args) != 5:    #4 arguments
        return helpMsg

    path = args[1]
    pkl_pref = args[2]
    no_partI = int(args[3])
    no_partIpO = int(args[4])

    ls_sc, ls_AF, ls_pos, ls_geno, ls_variant, ls_var_pos, ls_theta, ls_rho, ls_Ne = np.array([]), np.array([]), [], [], np.empty((0, CHRS)), [], np.array([]), np.array([]), np.array([])

    glob_thr = 0
    for ii in range(no_partI):
        with open(path+'/'+pkl_pref+'_'+str(ii)+'.pkl', 'rb') as f:
            list_sel_coef, list_freq, list_variant, _, list_geno, list_pos, list_variant_pos, _, theta, rho, Ne, _ = pickle.load(f)
            #print(type(list_sel_coef), type(list_freq), type(list_variant), type(list_geno), type(list_pos), type(list_variant_pos), type(theta), type(rho), type(Ne)); quit() # for debugging
            no_sims = len(list_sel_coef)
            ls_sc=np.append(ls_sc, list_sel_coef); ls_AF=np.append(ls_AF, list_freq); ls_pos+=list_pos; ls_geno+=list_geno; ls_variant=np.concatenate((ls_variant, list_variant)); ls_var_pos+=list_variant_pos; ls_theta=np.append(ls_theta, theta); ls_rho=np.append(ls_rho, rho); ls_Ne=np.append(ls_Ne, Ne)
            print("%% Buffering:", pkl_pref+'_'+str(ii)+'.pkl', no_sims, "to", glob_thr, flush=True)

        if (ii+1)%no_partIpO == 0:
            with open(pkl_pref+'/'+pkl_pref+'_'+str(glob_thr)+'.pkl', 'wb') as f:
                pickle.dump((ls_sc, ls_AF, ls_pos, ls_geno, ls_variant, ls_var_pos, ls_theta, ls_rho, ls_Ne), f, pickle.HIGHEST_PROTOCOL)
            print("%% Output to global thread", glob_thr, flush=True)
            print("%%", ls_sc.shape, ls_AF.shape, len(ls_pos), len(ls_geno), ls_variant.shape, len(ls_var_pos), ls_theta.shape, ls_rho.shape, ls_Ne.shape, flush=True)
            ls_sc, ls_AF, ls_pos, ls_geno, ls_variant, ls_var_pos, ls_theta, ls_rho, ls_Ne = np.array([]), np.array([]), [], [], np.empty((0, CHRS)), [], np.array([]), np.array([]), np.array([])
            glob_thr += 1

    return 0

sys.exit(main(sys.argv))
