#!/usr/bin/env python3

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "..")) 

from contextlib import contextmanager
import numpy as np
import pickle

import feature_util as fea

helpMsg = '''
        usage: $./sim2arg.py <no_st> <handle> <Ne> <thread #> <TAG>
            Takes a *partitioned* pickle file, runs RELATE to infer ARGs and extract features.
            - <no_st> >=0 : number of flanking gene trees to include on EACH side for feature extraction
            Thread # must match the partitioned pickle files
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

def main(args):
    if len(args) != 6:    #5 arguments
        return helpMsg

    no_st = int(args[1])
    handle = args[2]
    thread = args[4] # only used as a string
    tag = args[5]
    #mode = args[6]

    #Ne = args[3]
    Ne = int(args[3])
    Ne = str(2*Ne)

    pkl_path = handle+'/'+handle+'_pgv_'+thread+'.pkl'

    with open(pkl_path, 'rb') as f:  # Python 3: open(..., 'rb')
        list_pos, list_geno, list_variant, list_var_pos, iHS_df = pickle.load(f)

    lin_ls = []
    stat_ls = []
    relate_ls = [] # log-10 p-value for selection by RELATE

    wd = handle+'/'+handle+'_'+tag+'_temp_'+thread
    os.mkdir(wd, 0o755)

    with cd(wd):
        #iHS_idx = iHS_df[:, 0]-iHS_df[0, 0]
        for r_idx in range(len(list_pos)):
            #norm_iHS = iHS_df[iHS_idx==r_idx, 1:3]
            norm_iHS = None
            lin, stat, relate_p = fea.infer_ARG_fea(list_pos[r_idx], list_geno[r_idx], list_variant[r_idx], list_var_pos[r_idx], Ne, no_st, norm_iHS)
            lin_ls.append(lin)
            stat_ls.append(stat)
            relate_ls.append(relate_p)

    with open(handle+'/'+handle+'_'+tag+'_inf_fea_'+thread+'.pickle', 'wb') as f:
        pickle.dump((lin_ls, stat_ls, relate_ls), f, pickle.HIGHEST_PROTOCOL)

    os.rmdir(wd)

    return 0

sys.exit(main(sys.argv))
