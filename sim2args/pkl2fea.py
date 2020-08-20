#!/usr/bin/env python3

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "..")) 

from contextlib import contextmanager
import numpy as np
import pickle

from RELATE_util import run_RELATE
import utils

helpMsg = '''
        usage: $./pkl2fea.py <handle> <thread #> <TAG>
            Takes a *partitioned* pickle file, runs RELATE to infer ARGs and extract features.
            - <no_st> >=0 : number of flanking gene trees to include on EACH side for feature extraction
            Thread # must match the partitioned pickle files
'''

## Hyper-parameters ##

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
    if len(args) != 4:    #3 arguments
        return helpMsg

    handle = args[1]
    thread = args[2] # only used as a string
    tag = args[3]

    # Ne = int(args[3])
    # Nex2 = str(2*Ne) # commented out when Ne is part of pkl df
    # rho = float(args[4]) # deprecated

    pkl_path = handle+'/'+handle+'_'+thread+'.pkl'

    with open(pkl_path, 'rb') as f:  # Python 3: open(..., 'rb')
        ls_sc, ls_AF, list_pos, list_geno, list_variant, list_var_pos, list_theta, list_rho, list_Ne = pickle.load(f)

    tasks = len(ls_sc)

    success = np.empty(tasks, dtype=int)
    lin_df = np.empty((tasks, 2*no_ft+2, K), dtype=int)
    relate_ls = np.empty(tasks, dtype=float) # log-10 p-value for selection by RELATE

    wd = handle+'/'+handle+'_'+tag+'_temp_'+thread
    os.mkdir(wd, 0o755)

    with cd(wd):
        for r_idx in range(tasks):
            Ne = list_Ne[r_idx]
            Nex2 = str(2*Ne)
            mu = "{:.2e}".format(list_theta[r_idx]/(4*Ne*1e5)) # convert to string
            rho = list_rho[r_idx]/(4*Ne*1e5)*1e8 # in 1e8 unit
            lin_mtx, relate_p = inf_fea(list_pos[r_idx], list_geno[r_idx], list_variant[r_idx], list_var_pos[r_idx], Nex2, mu, rho)

            if lin_mtx is None:
                success[r_idx] = 0
            else:
                success[r_idx] = 1
                lin_df[r_idx] = lin_mtx
                relate_ls[r_idx] = relate_p

    np.savez_compressed(handle+'/'+handle+'_'+tag+'_inf_fea_'+thread, SC=ls_sc, CAF=ls_AF, fea_Mtx=lin_df, relateP=relate_ls, SUC=success)
    print("npz fmt - SC:CAF:fea_Mtx:relateP:SUC")
    print(f"Done_{np.sum(success)}/{success.shape}:", len(ls_sc), len(ls_AF), lin_df.shape, relate_ls.shape)

    os.rmdir(wd)

    return 0

sys.exit(main(sys.argv))
