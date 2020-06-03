#!/usr/bin/env python3

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from contextlib import contextmanager
import numpy as np
import feature_util as fea

## deprecated, use true tree for power analysis ##
helpMsg = '''
    usage: $./discoal2fea.py <`-s` or `-n`> <handle> <no_repl> <no_partition> <thr>

    Parse all discoal_<handle>/discoal_<handle>_*.discoal files
    Run RELATE and extract features
    Outputs
            discoal_<handle>/discoal_<handle>_inf_fea_*.npz

                SC, CAF, fea_Mtx

    <no_repl>: Total number of simulation runs
    <np_partition>: Number of partitions
    <thr>: Current thread #
    -s for sweep simulations
    -n for neutral simulations
'''

# @contextmanager is just an easier way of saying cd = contextmanager(cd)
@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def parse_hdswp_cmd(cmd_line):
    
    args = cmd_line.split()
    N_curr = int(args[args.index('-N') + 1])
    alpha = float(args[args.index('-a') + 1])
    AF = float(args[args.index('-c') + 1])
    selcoef = alpha/(2*N_curr)
    
    return selcoef, AF

def samp_var(geno_mtx, pos_ls, AF_min=0.2, AF_max=0.9):
    AF = np.mean(geno_mtx, axis=1)
    target = np.random.uniform(AF_min, AF_max)
    for idx in np.argsort(np.abs(AF - target)):
        if pos_ls[idx] > 0.4 and pos_ls[idx] < 0.6: # pick variant close to center (0.5)
            return pos_ls[idx], geno_mtx[idx], AF[idx]

def discoal_gt2fea(discoal_file, cat, no_ft):

    with open(discoal_file, "r") as discoalF:
        if cat == '-s':
            SC, CAF = parse_hdswp_cmd(discoalF.readline().strip())
        elif cat == '-n':
            SC = 0

        read_GT = False
        for line in discoalF:
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
            foc_var_pos, foc_var_gt, CAF = samp_var(gtm, var_pos, 0.01, 0.99)
        elif cat == '-s':
            foc_var_pos = 0.5
            foc_var_gt = gtm[var_pos == 0.5]

        feaMtx, _, _ = fea.infer_ARG_fea(var_pos, gtm, foc_var_gt, foc_var_pos, str(2*10000), "2.5e-8", 1.25, no_ft, None)
    
    return SC, CAF, feaMtx

def main(args):
    if len(args) != 6:    #5 argument
        return helpMsg
    
    mode = args[1]
    handle = args[2]
    no_sims = int(args[3])
    no_part = int(args[4])
    thr = int(args[5]) # zero-based
    
    ## Hyper-parameters ##
    K = 100 # number of time points for feature extraction (remember to change header of fea_util.py!)
    no_ft = 2 # number of flanking windows on each side for fea. extr.
    
    tasks = no_sims//no_part
    a_idx = thr*tasks
    if thr == no_part-1:
        b_idx = no_sims
    else:
        b_idx = (thr+1)*tasks

    SC_arr = np.empty(b_idx-a_idx) # of floats
    CAF_arr = np.empty(b_idx-a_idx) # of floats
    feaMtx_arr = np.empty((b_idx-a_idx, 2*no_ft+2, K)) # 3 dimensional

    wd = 'discoal_'+handle+'/discoal_'+handle+"_temp_"+str(thr)
    os.mkdir(wd, 0o755)

    with cd(wd):
        for samp in range(a_idx, b_idx):
            discoalF_path = '../discoal_'+handle+'_'+str(samp)+'.discoal'
            #if os.stat(discoalF_path).st_size == 0: continue
            SC, CAF, feaMtx = discoal_gt2fea(discoalF_path, mode, no_ft)

            SC_arr[samp-a_idx] = SC
            CAF_arr[samp-a_idx] = CAF

            feaMtx_arr[samp-a_idx, :, :] = feaMtx

    np.savez_compressed("discoal_"+handle+"/discoal_"+handle+"_inf_fea_"+str(thr), SC=SC_arr, CAF=CAF_arr, fea_Mtx=feaMtx_arr)
    print("Done:", SC_arr.shape, CAF_arr.shape, feaMtx_arr.shape)

    os.rmdir(wd)

sys.exit(main(sys.argv))
