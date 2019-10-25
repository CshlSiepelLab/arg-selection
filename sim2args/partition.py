#!/usr/bin/env python3

import sys, os
import pickle

helpMsg = '''
        usage: $./partition.py <pkl_path> <out_pref> <no_threads>
            - makes directory <out_pref>
            - outputs <out_pref>_pgv_<thread>.pickle files in directory
'''

def main(args):
    if len(args) != 4:    #3 arguments
        return helpMsg

    pkl_path = args[1]
    out_pref = args[2]
    no_threads = int(args[3])

    with open(pkl_path, 'rb') as f:  # Python 3: open(..., 'rb')
        list_sel_coef, list_freq, list_variant, list_tree, list_geno, list_pos = pickle.load(f)
        
    tasks = len(list_pos)//no_threads

    os.mkdir(out_pref, 0o755)

    for thread in range(1, no_threads+1):
        a_idx = (thread-1)*tasks
        if thread == no_threads:
            b_idx = len(list_pos)
        else:
            b_idx = thread*tasks
        with open(out_pref+'/'+out_pref+'_pgv_'+str(thread)+'.pickle', 'wb') as f:
            pickle.dump((list_pos[a_idx:b_idx], list_geno[a_idx:b_idx], list_variant[a_idx:b_idx]), f, pickle.HIGHEST_PROTOCOL)

    return 0

sys.exit(main(sys.argv))