#!/usr/bin/env python3

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "..")) 

from contextlib import contextmanager
import numpy as np
import pickle

from RELATE_util import run_RELATE
import utils

helpMsg = '''
        usage: $./pkl2fea.py <pkl_pref> <thread #> <tot_thr> <tot_files> <oTAG>
            Takes raw pickle files, runs RELATE to infer ARGs and extract features.
            - <no_st> >=0 : number of flanking gene trees to include on EACH side for feature extraction
            Smart threading: thread and file number must both be 0-indexed!!
            Make sure "tmp" and "out_threads" directories exist
'''

## Hyper-parameters ##

# bird - 24/54; CEU - 198; GBR - 182
CHRS=198 # number of chromosomes in simulation, needes for preallocation purposes

#discretT = np.logspace(1, 6, num=100)
time_file_path =  os.path.join(os.path.dirname(__file__), "..", "time.txt")
discretT = np.loadtxt(time_file_path)
discretT = discretT.astype(int)
K = len(discretT)
no_ft = 2 # of flanking gene trees to include on EACH side for feature extraction

# @contextmanager is just an easier way of saying cd = contextmanager(cd)
@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def inf_fea(pos_ls, geno_mtx, put_sel_var, var_pos, Ne, mu_string, rho_1e8):
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

    ts_inferred, RELATE_pval = run_RELATE(p, geno_mtx, Ne, var_ppos, rho_1e8, mu_string)
    if ts_inferred is None:
        return None, None

    trees = []
    intervals = np.empty((0,2), int)
    for tree in ts_inferred.trees():
        left, right = map(int, tree.interval)
        intervals = np.vstack((intervals, [left, right]))
        trees = np.append(trees, tree.newick(precision=5))

    intervals[-1, 1] = 1e5 # force last tree to cover the rest of the region

    lin_mtx = utils.ARG2feature(intervals, trees, var_ppos, put_sel_var, no_ft, discretT, 1)

    return lin_mtx, RELATE_pval # dim(lin_mtx)=(2*no_ft+2, K)

def main(args):
    if len(args) != 6:    #5 arguments
        return helpMsg

    handle = args[1]
    thread = int(args[2]) # xx only used as a string xx
    tot_thr = int(args[3])
    tot_files = int(args[4])
    tag = args[5]

    # Ne = int(args[3])
    # Nex2 = str(2*Ne) # commented out when Ne is part of pkl df
    # rho = float(args[4]) # deprecated

    if tot_thr < tot_files: # merge routine
        no_filePthr = tot_files//tot_thr
        file_idx0, file_idx1 = thread*no_filePthr, (thread+1)*no_filePthr
        print(f"++ MERGE ROUTINE THR{thread}: [{file_idx0}, {file_idx1})", flush=True)
        ls_sc, ls_AF, ls_pos, ls_geno, ls_variant, ls_var_pos, ls_theta, ls_rho, ls_Ne = np.array([]), np.array([]), [], [], np.empty((0, CHRS)), [], np.array([]), np.array([]), np.array([])
        
        for ii in range(file_idx0, file_idx1):
            with open(handle+'_'+str(ii)+'.pkl', 'rb') as f:
                list_sel_coef, list_freq, list_variant, _, list_geno, list_pos, list_variant_pos, _, theta, rho, Ne, _ = pickle.load(f)
                #print(type(list_sel_coef), type(list_freq), type(list_variant), type(list_geno), type(list_pos), type(list_variant_pos), type(theta), type(rho), type(Ne)); quit() # for debugging
            no_sims = len(list_sel_coef)
            ls_sc=np.append(ls_sc, list_sel_coef); ls_AF=np.append(ls_AF, list_freq); ls_pos+=list_pos; ls_geno+=list_geno; ls_variant=np.concatenate((ls_variant, list_variant)); ls_var_pos+=list_variant_pos; ls_theta=np.append(ls_theta, theta); ls_rho=np.append(ls_rho, rho); ls_Ne=np.append(ls_Ne, Ne)
            print("++ Buffering:", handle+'_'+str(ii)+'.pkl', no_sims, flush=True)

    else: # split routine
        no_thrPfile = tot_thr//tot_files
        file_idx = thread//no_thrPfile
        partition = thread%no_thrPfile

        print(f"++ SPLIT ROUTINE THR{thread}: pkl#{file_idx}, part#{partition}", flush=True)

        with open(handle+'_'+str(file_idx)+'.pkl', 'rb') as f:
            list_sel_coef, list_freq, list_variant, _, list_geno, list_pos, list_variant_pos, _, theta, rho, Ne, _ = pickle.load(f)
        # Partition
        no_sims = len(list_sel_coef)
        tasks = no_sims//no_thrPfile
        a_idx = partition*tasks
        if partition == no_thrPfile-1:
            b_idx = no_sims
        else:
            b_idx = (partition+1)*tasks
        print("++ Subsetting", handle+'_'+str(file_idx)+'.pkl', no_sims, ":", f"[{a_idx}, {b_idx})", flush=True)

        ls_sc, ls_AF, ls_pos, ls_geno, ls_variant, ls_var_pos, ls_theta, ls_rho, ls_Ne = list_sel_coef[a_idx:b_idx], list_freq[a_idx:b_idx], list_pos[a_idx:b_idx], list_geno[a_idx:b_idx], list_variant[a_idx:b_idx], list_variant_pos[a_idx:b_idx], theta[a_idx:b_idx], rho[a_idx:b_idx], Ne[a_idx:b_idx]

    print("++ Shapes:", len(ls_sc), len(ls_AF), len(ls_pos), len(ls_geno), ls_variant.shape, len(ls_var_pos), len(ls_theta), len(ls_rho), len(ls_Ne), flush=True)

    tasks = len(ls_sc)

    success = np.empty(tasks, dtype=int)
    lin_df = np.empty((tasks, 2*no_ft+2, K), dtype=int)
    relate_ls = np.empty(tasks, dtype=float) # log-10 p-value for selection by RELATE

    wd = 'tmp/'+tag+'_temp_'+str(thread)
    os.mkdir(wd, 0o755)

    with cd(wd):
        for r_idx in range(tasks):
            Ne = ls_Ne[r_idx]
            Nex2 = str(2*Ne)
            mu = "{:.2e}".format(ls_theta[r_idx]/(4*Ne*1e5)) # convert to string
            rho = ls_rho[r_idx]/(4*Ne*1e5)*1e8 # in 1e8 unit
            lin_mtx, relate_p = inf_fea(ls_pos[r_idx], ls_geno[r_idx], ls_variant[r_idx], ls_var_pos[r_idx], Nex2, mu, rho)

            if lin_mtx is None:
                success[r_idx] = 0
            else:
                success[r_idx] = 1
                lin_df[r_idx] = lin_mtx
                relate_ls[r_idx] = relate_p

    np.savez_compressed('out_threads/'+tag+'_inf_fea_'+str(thread), SC=ls_sc, CAF=ls_AF, fea_Mtx=lin_df, relateP=relate_ls, SUC=success)
    print("++ npz fmt - SC:CAF:fea_Mtx:relateP:SUC")
    print(f"++ Done_{np.sum(success)}/{success.shape}:", len(ls_sc), len(ls_AF), lin_df.shape, relate_ls.shape)

    os.rmdir(wd)

    return 0

sys.exit(main(sys.argv))
