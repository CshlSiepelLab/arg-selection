#!/usr/bin/env python3

import sys
import pickle
import numpy as np

helpMsg = '''
        usage: $./combine.py <handle> <swp_handle> <neu_handle> <out_TAG> <total threads>
            Combine discoal_<handle>/discoal_<swp_handle>_inf_fea_*.pickle
                    discoal_<handle>_neutral/discoal_<neu_handle>_inf_fea_*.pickle
            into inf_fea_<handle>_<out_TAG>.pkl
'''

def main(args):
    if len(args) != 6:    #5 arguments
        return helpMsg

    handle = args[1]
    swp_handle = args[2]
    neu_handle = args[3]
    tag = args[4]
    no_threads = int(args[5])

    feature_ls_neu = []
    feature_ls_swp = []

    for t in range(1, no_threads+1):
        with open('discoal_'+handle+'_neutral/'+'discoal_'+neu_handle+'_inf_fea_'+str(t)+'.pickle', 'rb') as f:
            feature_ls_neu+=pickle.load(f)
        with open('discoal_'+handle+'/'+'discoal_'+swp_handle+'_inf_fea_'+str(t)+'.pickle', 'rb') as f:
            feature_ls_swp+=pickle.load(f)
    
    with open('inf_fea_'+handle+'_'+tag+'.pkl', 'wb') as f:
        pickle.dump((feature_ls_neu, feature_ls_swp), f, pickle.HIGHEST_PROTOCOL)

    return 0

sys.exit(main(sys.argv))