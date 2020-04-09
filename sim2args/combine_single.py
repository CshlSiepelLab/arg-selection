#!/usr/bin/env python3

import sys
import pickle
import numpy as np
import os.path

helpMsg = '''
        usage: $./combine_single.py <handle> <TAG> <total threads> <oTAG>
            Combine <handle>/<handle>_<TAG>_inf_fea_*.pickle
            into <handle>_inf_fea_<oTAG>.pkl
            To load the resulting pickle file, run:
                with open('<handle>_inf_fea_<oTAG>.pkl', 'rb') as f:
                    lin_ls, stat_ls, relate_ls = pickle.load(f)
'''

def main(args):
    if len(args) != 5:    #4 arguments
        return helpMsg

    handle = args[1]
    tag = args[2]
    no_threads = int(args[3])
    oTAG = args[4]

    lin_ls = []
    stat_ls = []
    relate_ls = []

    for t in range(no_threads):
        Path = handle+'/'+handle+'_'+tag+'_inf_fea_'+str(t)+'.pickle'
        with open(Path, 'rb') as f:
            lin, stat, relate = pickle.load(f)
            print(Path+":", len(lin), len(stat), len(relate))
            lin_ls += lin
            stat_ls += stat
            relate_ls += relate
    
    print("Total:", len(lin_ls), len(stat_ls), len(relate_ls))
    print(lin_ls[np.random.randint(len(lin_ls))].shape,
        stat_ls[np.random.randint(len(stat_ls))].shape, sep='\n')

    with open(handle+'_inf_fea_'+oTAG+'.pkl', 'wb') as f:
        pickle.dump((lin_ls, stat_ls, relate_ls), f, pickle.HIGHEST_PROTOCOL)

    return 0

sys.exit(main(sys.argv))