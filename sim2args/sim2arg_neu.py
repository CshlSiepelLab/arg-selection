#!/usr/bin/env python3

import sys, os
from contextlib import contextmanager
import numpy as np
import pickle

import feature_util_neu as fea

helpMsg = '''
        usage: $./sim2arg_neu.py <fea> <handle> <Ne> <thread #, 1-based> <TAG>
            Takes a *partitioned* pickle file, runs RELATE to infer ARGs and extract features.
            <fea> determines the features extracted:
                1: single gene tree at site of interest
                3: incorporate the immediate surrounding trees (N/A yet!)
                5: incorporate the two immediate surrounding trees (N?A yet!)
'''

fea_xtract_dict = {
    '1': fea.infer_ARG_fea,
    '3': fea.infer_ARG_fea_3,
    '5': fea.infer_ARG_fea_5
}

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

    fx_mode = args[1]
    handle = args[2]
    thread = args[4] # only used as a string
    tag = args[5]

    #Ne = args[3]
    Ne = int(args[3])
    Ne = str(2*Ne)

    pkl_path = handle+'/'+handle+'_pgv_'+thread+'.pickle'

    with open(pkl_path, 'rb') as f:  # Python 3: open(..., 'rb')
        list_pos, list_geno, list_variant, list_var_pos = pickle.load(f)

    feature_ls = []

    wd = handle+'/'+handle+'_'+tag+'_temp_'+thread
    os.mkdir(wd, 0o755)

    with cd(wd):
    	for r_idx in range(len(list_pos)):
    		feature_ls.append(fea_xtract_dict[fx_mode](list_pos[r_idx], list_geno[r_idx], list_variant[r_idx], list_var_pos[r_idx], Ne))

    with open(handle+'/'+handle+'_'+tag+'_inf_fea_'+thread+'.pickle', 'wb') as f:
    	pickle.dump(feature_ls, f, pickle.HIGHEST_PROTOCOL)

    os.rmdir(wd)

    return 0

sys.exit(main(sys.argv))
