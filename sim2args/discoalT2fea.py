#!/usr/bin/env python3

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

import dendropy
import numpy as np
import utils

helpMsg = '''
    usage: $./discoalT2fea.py <`-s` or `-n`> <handle> <no_repl> <no_partition> <thr> 

    Parse all discoal_<handle>/discoal_<handle>_*.discoal files
    Outputs discoal_<handle>/discoal_<handle>_tru_fea.npz (feature from true gene trees)

                SC, CAF, fea_Mtx

    <no_repl>: Total number of simulation runs
    <np_partition>: Number of partitions
    <thr>: Current thread #
    -s for sweep simulations
    -n for neutral simulations
'''

## Hyper-parameters ##
discreT = utils.discretizeT(0.005, 2e4, 100, 2*1e4)
K = len(discreT) # number of time points for feature extraction
no_ft = 2 # number of flanking windows on each side for fea. extr.

def parse_tree_line(line):

    intv, nwk = line.split(']')
    intv = int(intv[1:])
    #tr = dendropy.Tree.get(data=nwk, schema="newick")
    
    return intv, nwk

def parse_hdswp_cmd(cmd_line):
    
    args = cmd_line.split()
    N_curr = int(args[args.index('-N') + 1])
    alpha = float(args[args.index('-a') + 1])
    AF = float(args[args.index('-c') + 1])
    selcoef = alpha/(2*N_curr)
    
    return selcoef, AF

def process_discoalT(discoal_file, cat):

    # Pass: Read sc, var_pos, gt_mtx and trees
    with open(discoal_file, "r") as discoalF:
        if cat == '-s':
            SC, CAF = parse_hdswp_cmd(discoalF.readline().strip())
        elif cat == '-n':
            SC = 0

        read_GT = False
        intvl_ls = []
        nwk_ls = []
        for line in discoalF:
            if line[0] == '[':
                intvl, nwk = parse_tree_line(line.strip())
                intvl_ls.append(intvl)
                nwk_ls.append(nwk)
                continue
            if line[:8] == "segsites":
                segsites = int(line.strip().split()[1])
                gtm = np.empty((0, segsites), dtype=np.int8)
                continue
            if line[:9] == "positions":
                var_pos = np.array(line.strip().split()[1:], dtype=float)
                read_GT = True
                continue
            if read_GT:
                gtm = np.vstack((gtm, np.fromstring(line.strip(), dtype=np.int8) - ord("0")))

    gtm = np.transpose(gtm)
    if cat == '-n':
        samp_idx = utils.samp_var(gtm, var_pos, 0.01, 0.99, 0.4, 0.6)
        foc_var_pos = var_pos[samp_idx]
        foc_var_gt = gtm[samp_idx]
        CAF = np.mean(foc_var_gt)

    elif cat == '-s':
        foc_var_pos = 0.5
        foc_var_gt = gtm[var_pos == 0.5].flatten()

    intvl_end = np.cumsum(intvl_ls)
    intvl_df = np.stack((np.concatenate(([0], intvl_end[:-1])), intvl_end), axis=-1)
    # Feature extraction
    feaMtx = utils.ARG2feature(intvl_df, nwk_ls, int(foc_var_pos*1e5), foc_var_gt, no_ft, discreT, 0)
    
    return SC, CAF, feaMtx

def main(args):
    if len(args) != 6:    #5 argument
        return helpMsg
    
    mode = args[1]
    handle = args[2]
    no_sims = int(args[3])
    no_part = int(args[4])
    thr = int(args[5]) # zero-based

    tasks = no_sims//no_part
    a_idx = thr*tasks
    if thr == no_part-1:
        b_idx = no_sims
    else:
        b_idx = (thr+1)*tasks

    SC_arr = np.empty(b_idx-a_idx) # of floats
    CAF_arr = np.empty(b_idx-a_idx) # of floats
    feaMtx_arr = np.empty((b_idx-a_idx, 2*no_ft+2, K)) # 3 dimensional

    for samp in range(a_idx, b_idx):
        discoalF_path = 'discoal_'+handle+'/discoal_'+handle+'_'+str(samp)+'.discoal'
        #if os.stat(discoalF_path).st_size == 0: continue
        SC, CAF, feaMtx = process_discoalT(discoalF_path, mode)
        print(samp, "_Done_", flush=True)
        SC_arr[samp-a_idx] = SC
        CAF_arr[samp-a_idx] = CAF

        feaMtx_arr[samp-a_idx, :, :] = feaMtx

    np.savez_compressed("discoal_"+handle+"/discoal_"+handle+"_tru_fea_"+str(thr), SC=SC_arr, CAF=CAF_arr, fea_Mtx=feaMtx_arr)
    print("npz fmt - SC:CAF:fea_Mtx")
    print("Done:", SC_arr.shape, CAF_arr.shape, feaMtx_arr.shape)

sys.exit(main(sys.argv))
