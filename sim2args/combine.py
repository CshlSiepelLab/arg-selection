#!/usr/bin/env python3

import sys
import pickle
import numpy as np

helpMsg = '''
        usage: $./combine.py <swp_handle> <neu_handle> <TAG> <total threads>
            Combine <swp_handle>/<swp_handle>_<TAG>_inf_fea_*.pickle
                    <neu_handle>/<neu_handle>_<TAG>_inf_fea_*.pickle
            into inf_fea_<swp_handle>_<TAG>.pkl
            To load the resulting pickle file, run:
                with open('inf_fea_<swp_handle>_<TAG>.pkl', 'rb') as f:
                    lin_ls_neu, stat_ls_neu, lin_ls_swp, stat_ls_swp = pickle.load(f)
'''

def main(args):
    if len(args) != 5:    #4 arguments
        return helpMsg

    #handle = args[1]
    swp_handle = args[1]
    neu_handle = args[2]
    tag = args[3]
    no_threads = int(args[4])

    lin_ls_neu = []
    stat_ls_neu = []
    
    lin_ls_swp = []
    stat_ls_swp = []

    for t in range(1, no_threads+1):
        with open(neu_handle+'/'+neu_handle+'_'+tag+'_inf_fea_'+str(t)+'.pickle', 'rb') as f:
            lin, stat = pickle.load(f)
            lin_ls_neu += lin
            stat_ls_neu += stat

        with open(swp_handle+'/'+swp_handle+'_'+tag+'_inf_fea_'+str(t)+'.pickle', 'rb') as f:
            lin, stat = pickle.load(f)
            lin_ls_swp += lin
            stat_ls_swp += stat
    
    print(len(lin_ls_neu), len(stat_ls_neu), len(lin_ls_swp), len(stat_ls_swp))
    print(lin_ls_neu[np.random.randint(len(lin_ls_neu))].shape,
        stat_ls_neu[np.random.randint(len(stat_ls_neu))].shape,
        lin_ls_swp[np.random.randint(len(lin_ls_swp))].shape,
        stat_ls_swp[np.random.randint(len(stat_ls_swp))].shape, sep='\n')

    with open('inf_fea_'+swp_handle+'_'+tag+'.pkl', 'wb') as f:
        pickle.dump((lin_ls_neu, stat_ls_neu, lin_ls_swp, stat_ls_swp), f, pickle.HIGHEST_PROTOCOL)

    return 0

sys.exit(main(sys.argv))