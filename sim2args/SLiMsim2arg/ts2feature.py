#!/usr/bin/env python3

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../..")) # directory of module utils

import numpy as np
import dendropy
import tskit
import tszip

import utils

# MACRO: simulated region [0, 1e5)
win_l = np.array([-4.5e4, -3.5e4, -2.5e4, -1.5e4, -1e4, -5e3, 0, 1, 1+5e3, 1+1e4, 1+1.5e4, 1+2.5e4, 1+3.5e4])
win_r = np.array([-3.5e4, -2.5e4, -1.5e4, -1e4, -5e3, 0, 1, 1+5e3, 1+1e4, 1+1.5e4, 1+2.5e4, 1+3.5e4, 1+4.5e4])
discreT = utils.discretizeT()
K = len(discreT)
W = len(win_l)

helpMsg = '''
        usage: $./ts2feature.py <`n`/`s`> <no_sims/meta_file_path> <thr> <tot_thr> <in_pref> <tree_type> <outPref>
        meta_file is .npz file
        <in_pref> should be the complete prefix, including directory and file name prefix
        <tree_type> example: `tru.trees`, `inf.trees` or `inf.trees.tsz`
        (e.g. SLiM_trial_neu_ts/SLiM_trial_neu_ts)
'''

def main(args):
    if len(args) != 8:    #7 arguments
        return helpMsg

    mode = args[1]
    if mode == 'n':
        no_sims = int(args[2])
    elif mode == 's':
        metaF = args[2]
    else:
        return helpMsg

    thr = int(args[3]) # should be 1-indexed
    tot_thr = int(args[4])
    inPref = args[5]
    tree_type = args[6]
    outPref = args[7]

    if mode == 's':
        with np.load(metaF) as meta:
            idx_ls = meta["idx"]
            sc_ls = meta["sc"]
            onset_ls = meta["onset"]
            caf_ls = meta["caf"]
        no_sims = idx_ls.shape[0]

    tasks = no_sims//tot_thr
    a_idx = (thr-1)*tasks # inclusive
    if thr == tot_thr:
        b_idx = no_sims # exclusive
    else:
        b_idx = thr*tasks # exclusive
    # indices are 0-based

    print("Processing: [", a_idx, b_idx, ")", flush=True)

    fea_df = np.empty((0, W*3, K)) # dataframe containing the features
    meta_df = np.empty((0, 5)) # dataframe containing the meta-data [idx, sc, onset, AF, var_pos]

    for r_idx in range(a_idx, b_idx):
        if mode == 'n':
            ID = r_idx + 1 # convert 0-based index to 1-based index
            ts_path = inPref+"_"+str(ID)+"_"+tree_type
            if not os.path.isfile(ts_path): continue

            if tree_type[-3:] == "tsz":
                ts_eg = tszip.decompress(ts_path)
            else:
                ts_eg = tskit.load(ts_path)

            GTM = ts_eg.genotype_matrix()
            var_pos = utils.get_site_ppos(ts_eg)

            samp_idx = utils.samp_var(GTM, var_pos)
            vOI_gt = GTM[samp_idx, :].flatten()
            vOI_pos = var_pos[samp_idx]
            sc = 0
            onset = 0
            AF = np.mean(vOI_gt)
            
        elif mode == 's':
            ID = idx_ls[r_idx] # retrieve 1-based index from metadata
            ts_path = inPref+"_"+str(ID)+"_"+tree_type

            if tree_type[-3:] == "tsz":
                ts_eg = tszip.decompress(ts_path)
            else:
                ts_eg = tskit.load(ts_path)

            GTM = ts_eg.genotype_matrix()
            var_pos = utils.get_site_ppos(ts_eg)

            vOI_gt = GTM[var_pos==50000, :].flatten()
            vOI_pos = 50000
            sc = sc_ls[r_idx]
            onset = onset_ls[r_idx]
            AF = caf_ls[r_idx]

        print(">>>", r_idx, "/", b_idx, ts_path, ":", vOI_pos, flush=True)
        fea_mtx = utils.ts2feature(ts_eg, vOI_pos, vOI_gt, win_l, win_r, discreT)

        fea_df = np.concatenate((fea_df, np.array([fea_mtx])))
        meta_df = np.vstack((meta_df, [ID, sc, onset, AF, vOI_pos]))

    print(fea_df.shape, meta_df.shape)
    np.save(outPref+"_fea_"+str(thr), fea_df)
    np.save(outPref+"_meta_"+str(thr), meta_df)

    return 0

sys.exit(main(sys.argv))
