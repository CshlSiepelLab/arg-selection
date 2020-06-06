#!/usr/bin/env python3

import sys
import numpy as np
import os.path

helpMsg = '''
        usage: $./combine_fea_npz.py <prefix> <KEY> <base> <total threads> <oTAG>
            <base> (0- zero-base; 1- one-base)
            Combine <prefix>_*.npz ==> <oTAG>.npz
            <KEY> describes the key to arrays in the .npz files
                e.g. SC:CAF:fea_Mtx (colon separated)
'''

def main(args):
    if len(args) != 6:    #5 arguments
        return helpMsg

    pref = args[1]
    keys = args[2].split(":")
    base = int(args[3])
    no_threads = int(args[4])
    oTAG = args[5]

    no_dfs = len(keys)
    dfs = {}

    for thr in range(no_threads):
        t = thr + base
        Path = pref+'_'+str(t)+'.npz'
        if os.path.isfile(Path):
            with np.load(Path) as data:
                for key in keys:
                    if t == base:
                        dfs[key] = data[key]
                    else:
                        dfs[key] = np.concatenate((dfs[key], data[key]))
        else:
            print(Path, ": file does not exist, SKIPPING")
    
    np.savez_compressed(oTAG, **dfs)
    ## Print meta-data ##
    for key in keys:
        print(f"KEY:`{key}`\t\tSHAPE: {dfs[key].shape}")

    return 0

sys.exit(main(sys.argv))