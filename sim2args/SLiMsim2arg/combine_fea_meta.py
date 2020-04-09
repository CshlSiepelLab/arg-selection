#!/usr/bin/env python3

import sys
import numpy as np
import os.path

helpMsg = '''
        usage: $./combine_fea_meta.py <PATH_pref> <total threads> <oTAG>
            Combine <PATH_pref>_meta_*.npy => <oTAG>_meta.npy
                    <PATH_pref>_fea_*.npy  => <oTAG>_fea.npy
'''

## MANUALLY DEFINE MACRO ##
K = 100
W = 13

def main(args):
    if len(args) != 4:    #3 arguments
        return helpMsg

    path_pref = args[1]
    no_threads = int(args[2])
    oTAG = args[3]

    fea_df = np.empty((0, W*3, K)) # dataframe containing the features
    meta_df = np.empty((0, 5)) # dataframe containing the meta-data [idx, sc, onset, AF, var_pos]

    for t in range(1, no_threads+1):
        fea_df = np.concatenate((fea_df, np.load(path_pref+'_fea_'+str(t)+'.npy')))
        meta_df = np.concatenate((meta_df, np.load(path_pref+'_meta_'+str(t)+'.npy')))

    np.save(oTAG+'_fea', fea_df)
    np.save(oTAG+'_meta', meta_df)

    print(path_pref+':', fea_df.shape, meta_df.shape)

    return 0

sys.exit(main(sys.argv))