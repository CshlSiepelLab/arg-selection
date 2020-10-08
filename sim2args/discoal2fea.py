#!/usr/bin/env python3

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from contextlib import contextmanager
import numpy as np
import pickle
from RELATE_util import run_RELATE
import utils

helpMsg = '''
    usage: $./discoal2fea.py <handle> <thr> <no_pkl> <no_thr>

    Parse all discoal_pkl/discoal_<handle>_*.pkl files
    Run RELATE and extract features
    Outputs
            discoal_inf_fea/discoal_<handle>_inf_fea_*.npz

                SC_arr, SAF_arr, CAF_arr, onset_arr, fea_Mtx_df

    <thr>: Current thread #
    -s for sweep simulations
    -n for neutral simulations
'''

# MACRO: simulated region [0, 1e5)
win_l = np.array([-4.5e4, -3.5e4, -2.5e4, -1.5e4, -1e4, -5e3, 0, 1, 1+5e3, 1+1e4, 1+1.5e4, 1+2.5e4, 1+3.5e4])
win_r = np.array([-3.5e4, -2.5e4, -1.5e4, -1e4, -5e3, 0, 1, 1+5e3, 1+1e4, 1+1.5e4, 1+2.5e4, 1+3.5e4, 1+4.5e4])
discreT = utils.discretizeT()
K = len(discreT)
W = len(win_l)

# @contextmanager is just an easier way of saying cd = contextmanager(cd)
@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def inf_fea(put_sel_var_pos, pos_ls, geno_mtx, Nex2=str(2*1e8), mu_string="3e-8", rho_1e8=0.7):
    '''Format input, run RELATE on variants of a simulated region and extract features of a region from inferred tree sequence
    mu and rho are rounded for both teo and LR
    '''

    p = utils.discrete_pos(pos_ls)

    ts_inferred, _ = run_RELATE(p, geno_mtx, Nex2, -1, rho_1e8, mu_string) # ignore relate selection inference
    if ts_inferred is None:
        return None

    foc_var_gt = geno_mtx[pos_ls == put_sel_var_pos].flatten()
    foc_var_ppos = p[pos_ls == put_sel_var_pos][0] # BUG FIX: convert decimal position to int!!
    lin_mtx = utils.ts2feature(ts_inferred, foc_var_ppos, foc_var_gt, win_l, win_r, discreT)

    return lin_mtx

def main(args):
    if len(args) != 5:    #4 argument(s)
        return helpMsg
    
    #mode = args[1]
    handle = args[1]
    thr = int(args[2]) # one-based
    no_pkl = int(args[3]) # pkl idx is one-based
    no_thr = int(args[4]) # should be multiples of no_pkl

    thr_pp = no_thr//no_pkl
    pkl_idx = (thr-1)//thr_pp+1

    with open(f"discoal_pkl/discoal_{handle}_{pkl_idx}.pkl", "rb") as f:
        SC_arr, SAF_arr, CAF_arr, onset_arr, foc_var_pos_arr, var_pos_ls, gtm_ls = pickle.load(f)

    no_pkl_sims = SC_arr.shape[0]
    sims_p_thr = no_pkl_sims//thr_pp
    a_idx = ((thr-1)%thr_pp)*sims_p_thr
    if thr%thr_pp == 0:
        b_idx = no_pkl_sims
    else:
        b_idx = a_idx + sims_p_thr
    no_sims = b_idx - a_idx
    feaMtx_df = np.empty((no_sims, W*3, K)) # dataframe containing the features
    success = np.empty(no_sims, dtype=int)

    print(f"%Thread {thr}, pkl file {pkl_idx}:{a_idx}-{b_idx}", flush=True)

    wd = 'discoal_temp/discoal_'+handle+"_temp_"+str(thr)
    os.mkdir(wd, 0o755)

    with cd(wd):
        for samp in range(a_idx, b_idx):
            #if os.stat(discoalF_path).st_size == 0: continue
            feaMtx = inf_fea(foc_var_pos_arr[samp], var_pos_ls[samp], gtm_ls[samp])
            
            if feaMtx is None:
                success[samp-a_idx] = 0
                print(samp, "_FAILED_", flush=True)
            else:
                success[samp-a_idx] = 1
                feaMtx_df[samp-a_idx, :, :] = feaMtx
                print(samp, "_SUCCEEDED_", flush=True)

    np.savez_compressed("discoal_inf_fea/discoal_"+handle+"_inf_fea_"+str(thr),
        SC=SC_arr[a_idx:b_idx], SAF=SAF_arr[a_idx:b_idx], CAF=CAF_arr[a_idx:b_idx], onset=onset_arr[a_idx:b_idx], fea_Mtx=feaMtx_df, SUC=success)
    print("%npz fmt - SC:SAF:CAF:onset:fea_Mtx:SUC")
    print(f"%Done_{np.sum(success)}/{success.shape}:", feaMtx_df.shape)

    os.rmdir(wd)

sys.exit(main(sys.argv))
