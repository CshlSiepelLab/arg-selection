#!/usr/bin/env python3

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from contextlib import contextmanager
import numpy as np

from RELATE_util import run_RELATE
import utils

## Feature from inferred tree for power analysis ##
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

## Hyper-parameters ##
discreT = utils.discretizeT(0.005, 2e4, 100, 2*1e4)
K = len(discreT) # number of time points for feature extraction
no_ft = 2 # number of flanking windows on each side for fea. extr.

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

def inf_fea(pos_ls, geno_mtx, put_sel_var, var_pos, Nex2, mu_string, rho_1e8):
    '''Format input, run RELATE on variants of a simulated region and extract features of a region from inferred tree sequence

        For neutral regions, `var_pos` is the site chosen to have a matching allele frequency;
        for sweep regions, `var_pos` is specified to be 0.5
    '''

    p = utils.discrete_pos(pos_ls)

    if var_pos == 0.5:
        var_ppos = 50000
    else:
        var_idx = np.where(pos_ls == var_pos)[0][0]
        var_ppos = p[var_idx]

    ts_inferred, _ = run_RELATE(p, geno_mtx, Nex2, -1, rho_1e8, mu_string) # ignore relate selection inference
    if ts_inferred is None:
        return None

    trees = []
    intervals = np.empty((0,2), int)
    for tree in ts_inferred.trees():
        left, right = map(int, tree.interval)
        intervals = np.vstack((intervals, [left, right]))
        trees = np.append(trees, tree.newick(precision=5))

    intervals[-1, 1] = 1e5 # force last tree to cover the rest of the region

    lin_mtx = utils.ARG2feature(intervals, trees, var_ppos, put_sel_var, no_ft, discretT, 1)

    return lin_mtx # dim(lin_mtx)=(2*no_ft+2, K)

def discoal_gt2fea(discoal_file, cat):

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
        samp_idx = utils.samp_var(gtm, var_pos, 0.01, 0.99, 0.4, 0.6)
        foc_var_pos = var_pos[samp_idx]
        foc_var_gt = gtm[samp_idx]
        CAF = np.mean(foc_var_gt)

    elif cat == '-s':
        foc_var_pos = 0.5
        foc_var_gt = gtm[var_pos == 0.5].flatten()

    feaMtx = inf_fea(var_pos, gtm, foc_var_gt, foc_var_pos, str(2*10000), "2.5e-8", 1.25)
    
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
    success = np.empty(b_idx-a_idx, dtype=int)

    wd = 'discoal_'+handle+'/discoal_'+handle+"_temp_"+str(thr)
    os.mkdir(wd, 0o755)

    with cd(wd):
        for samp in range(a_idx, b_idx):
            discoalF_path = '../discoal_'+handle+'_'+str(samp)+'.discoal'
            #if os.stat(discoalF_path).st_size == 0: continue
            SC, CAF, feaMtx = discoal_gt2fea(discoalF_path, mode)
            
            SC_arr[samp-a_idx] = SC
            CAF_arr[samp-a_idx] = CAF
            if feaMtx is None:
                success[samp-a_idx] = 0
                print(samp, "_FAILED_", flush=True)
            else:
                success[samp-a_idx] = 1
                feaMtx_arr[samp-a_idx, :, :] = feaMtx
                print(samp, "_SUCCEEDED_", flush=True)

    np.savez_compressed("discoal_"+handle+"/discoal_"+handle+"_inf_fea_"+str(thr), SC=SC_arr, CAF=CAF_arr, fea_Mtx=feaMtx_arr, SUC=success)
    print("npz fmt - SC:CAF:fea_Mtx:SUC")
    print(f"Done_{np.sum(success)}/{success.shape}:", SC_arr.shape, CAF_arr.shape, feaMtx_arr.shape)

    os.rmdir(wd)

sys.exit(main(sys.argv))
