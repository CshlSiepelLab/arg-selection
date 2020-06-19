#!/usr/bin/env python3

import sys, os
import numpy as np

helpMsg = '''
        usage: $./meta2npz.py <meta_file_path> <outPref>

        4th column of meta file must contain the directory prefix of simulation files
        RELATIVE TO CWD()
'''

def main(args):
    if len(args) != 3:    #2 arguments
        return helpMsg

    metaF = args[1]
    outPref = args[2]

    meta_data = np.genfromtxt(metaF, usecols=(1, 2, 3, 4), dtype=None)
    # e.g. `%% 0.100992 2347 0.520421 SLiM_trial_swp/SLiM_trial_swp_4501_temp` <= artifact of "temp"
    no_sims = meta_data.shape[0]

    idx_ls = []
    sc_ls = []
    onset_ls = []
    caf_ls = []

    cnt = 0
    for r_idx in range(no_sims):
        ID = int(meta_data[r_idx][3].split(b'_')[-2]) # retrieve 1-based index from meta file
        sim_path = meta_data[r_idx][3].decode()[:-4]+"samp.trees"

        if not os.path.isfile(sim_path): continue
        idx_ls.append(ID)
        sc_ls.append(meta_data[r_idx][0])
        onset_ls.append(meta_data[r_idx][1])
        caf_ls.append(meta_data[r_idx][2])
        cnt += 1

    np.savez_compressed(outPref+"_meta", idx=idx_ls, sc=sc_ls, onset=onset_ls, caf=caf_ls)
    print(cnt, "simulations exist")

    return 0

sys.exit(main(sys.argv))
