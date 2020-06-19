#!/usr/bin/env python3

import os, sys
import subprocess
sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
 
from contextlib import contextmanager
import numpy as np

from RELATE_util import run_RELATE
import utils

## Feature from inferred tree for power analysis ##
helpMsg = '''
    usage: $./discoal_sim2fea.py <`swp` or `neu`> <META> <NOCHR> <thr>

    Simulation de novo
    Run RELATE and extract features
    Outputs
            "discoal_inf_fea/discoal_inf_fea_{handle}_{no_chrs}_{mode}_{thr}.npz"

                SC, CAF, fea_Mtx

'''

DISCOAL_PATH = '/sonas-hs/siepel/hpc_norepl/home/mo/discoal-master'

## Simulation parameters ##
Ne = 10000 # constant population size (10000)
pos = 0.5 #beneficial mutation positions (0-1)
length = 100000 #length of the simulated region
mu = 2.5e-8
rho = 1.25e-8
SIGMA = 40 # time discretization for sweep simulation range: (4, 400),dt=(1/sigma*N)

theta = 4*Ne*mu*length
R = 4*Ne*rho*length

## Feature parameters ##
discreT = utils.discretizeT(0.005, 2e4, 100) # RELATE tree branch langth in GENERATIONS!!
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

    lin_mtx = utils.ARG2feature(intervals, trees, var_ppos, put_sel_var, no_ft, discreT, 1)

    return lin_mtx # dim(lin_mtx)=(2*no_ft+2, K)

def discoal_gt2fea(discoal_file, cat):

    with open(discoal_file, "r") as discoalF:
        if cat == 'swp':
            SC, CAF = parse_hdswp_cmd(discoalF.readline().strip())
        elif cat == 'neu':
            SC = 0

        read_GT = False
        seek_onset = False
        onset_gen = -1
        for line in discoalF:
            if line[:8] == "segsites":
                segsites = int(line.strip().split()[1])
                gtm = np.empty((0, segsites), dtype=np.int8)
                continue
            if line[:9] == "positions":
                var_pos = np.array(line.strip().split()[1:], dtype=float)
                read_GT = True
                continue
            if line[:4] == "Freq":
                seek_onset = True
                continue
            if seek_onset:
                gbp_der_anc = line.strip().split()
                if len(gbp_der_anc) == 3:
                    onset_gen = float(gbp_der_anc[0]) # in coalc. unit
                elif onset_gen != -1:
                    seek_onset = False
                continue
            if read_GT:
                gtm = np.vstack((gtm, np.fromstring(line.strip(), dtype=np.int8) - ord("0")))

    gtm = np.transpose(gtm)
    if cat == 'neu':
        samp_idx = utils.samp_var(gtm, var_pos, 0.01, 0.99, 0.4, 0.6)
        foc_var_pos = var_pos[samp_idx]
        foc_var_gt = gtm[samp_idx]
        CAF = np.mean(foc_var_gt)

    elif cat == 'swp':
        foc_var_pos = 0.5
        foc_var_gt = gtm[var_pos == 0.5].flatten()

    feaMtx = inf_fea(var_pos, gtm, foc_var_gt, foc_var_pos, str(2*10000), "2.5e-8", 1.25)
    
    return SC, CAF, onset_gen, feaMtx

def main(args):
    if len(args) != 5:    #4 argument
        return helpMsg
    
    mode = args[1] # <`swp` or `neu`>
    handle = args[2]
    no_chrs = int(args[3])
    thr = int(args[4]) # one-based
    
    meta_df = np.genfromtxt(f"ID_SC_AF/discoal_{handle}_{thr}_id_sc_af.txt")
    no_sims = meta_df.shape[0]
    print(f"discoal_{handle}_{no_chrs}_{thr}: {no_sims} sims", flush=True)

    SC_arr = meta_df[:, 1] # of floats
    CAF_arr = meta_df[:, 2] # of floats
    onset_arr = np.empty(no_sims)
    feaMtx_arr = np.empty((no_sims, 2*no_ft+2, K)) # 3 dimensional
    success = np.empty(no_sims, dtype=int)

    wd = f"discoal_temp/discoal_{handle}_{no_chrs}_{mode}_{thr}"
    os.mkdir(wd, 0o755)
    temp_discoalF = f"discoal_{handle}_{no_chrs}_{mode}_{thr}_temp.discoal"

    with cd(wd):
        for samp in range(no_sims):
            #if os.stat(discoalF_path).st_size == 0: continue
            if mode == 'swp':
                sel = 2*Ne*SC_arr[samp]
                discoal_cmd = [f"{DISCOAL_PATH}/discoal", str(no_chrs), "1", str(length),
                "-t", str(theta), "-r", str(R), "-N", str(Ne),
                "-i", str(SIGMA), "-c", str(CAF_arr[samp]), "-ws", "0", "-a", str(sel), "-x", str(pos), "-T"]

            elif mode == 'neu':
                discoal_cmd = [f"{DISCOAL_PATH}/discoal", str(no_chrs), "1", str(length),
                "-t", str(theta), "-r", str(R), "-T"]

            loop_cnt = 0
            while True:
                loop_cnt += 1
                print(f"SIM:{samp}, ATTEMPT:{loop_cnt}", flush=True)
                with open(temp_discoalF, "w") as outF:
                    discoal_proc = subprocess.run(discoal_cmd, stdout=outF)
                if discoal_proc.returncode == 0: break

            SC, CAF, onset_gen, feaMtx = discoal_gt2fea(temp_discoalF, mode)
            
            if mode == 'neu':
                SC_arr[samp] = SC
                CAF_arr[samp] = CAF

            print(SC, CAF, ";", SC_arr[samp], CAF_arr[samp], ";", onset_gen) # sanity check
            onset_arr[samp] = onset_gen

            if feaMtx is None:
                success[samp] = 0
                print(samp, "_FAILED_", flush=True)
            else:
                success[samp] = 1
                feaMtx_arr[samp, :, :] = feaMtx
                print(samp, "_SUCCEEDED_", flush=True)

            os.remove(temp_discoalF)

    np.savez_compressed(f"discoal_inf_fea/discoal_inf_fea_{handle}_{no_chrs}_{mode}_{thr}", SC=SC_arr, CAF=CAF_arr, onset=onset_arr, fea_Mtx=feaMtx_arr, SUC=success)
    print("npz fmt - SC:CAF:onset:fea_Mtx:SUC")
    print(f"Done_{np.sum(success)}/{success.shape}:", SC_arr.shape, CAF_arr.shape, onset_arr.shape, feaMtx_arr.shape)

    os.rmdir(wd)

sys.exit(main(sys.argv))
