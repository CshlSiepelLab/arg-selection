#!/usr/bin/env python3

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "..")) 

import dendropy
import tskit
import pickle
import gzip
import numpy as np

import utils

helpMsg = '''
        usage: $./genFeatures_win.py <.trees PATH> <.meta PATH> <out prefix> <min_DAF>
'''

win_l = np.array([-4.5e4, -3.5e4, -2.5e4, -1.5e4, -1e4, -5e3, 0, 1, 1+5e3, 1+1e4, 1+1.5e4, 1+2.5e4, 1+3.5e4])
win_r = np.array([-3.5e4, -2.5e4, -1.5e4, -1e4, -5e3, 0, 1, 1+5e3, 1+1e4, 1+1.5e4, 1+2.5e4, 1+3.5e4, 1+4.5e4])
discreT = utils.discretizeT()
K = len(discreT)
W = len(win_l)

def main(args):
    if len(args) != 5:    #4 arguments
        return helpMsg

    ts_file=args[1]
    meta_file=args[2]
    out_pref=args[3]
    min_DAF=float(args[4])
    
    ts = tskit.load(ts_file)
    # meta = np.genfromtxt(meta_file, dtype=None, encoding=None) 
    # rs_ID | pos | code e.g. ('rs59993883', 135500045, 0)
    meta_ID = np.loadtxt(meta_file, dtype=str, usecols=0, encoding=None)
    meta_pos = np.loadtxt(meta_file, dtype=int, usecols=1)
    meta_code = np.loadtxt(meta_file, dtype=int, usecols=2)
    print(f"<> No. sites in: ts_file-{ts.num_sites}, meta_file-{meta_ID.shape[0]}; {meta_pos.shape[0]}; {meta_code.shape[0]}", flush=True)
    # ^^ so far not using meta file

    no_haptypes = ts.num_samples
    min_DAC = round(min_DAF*no_haptypes)
    print(f"<> Min DAC threshold: {min_DAC}/{no_haptypes}")
    
    # feature matrices for variant sites
    #id_ls = []
    daf_ls = []
    #gt_ls = []
    pos_ls = []
    lin_df = np.empty((0, W*3, K))
    #tree_ls = [] #nwk strings of the trees at each site

    # initialize
    intervals = np.empty((0, 2), int) # of row indefinite, total coverage must be 100k
    #nwk_tr_ls = []  # list of newick strings, # of row matches above

    leading_tr = ts.first()
    left, right = map(int, leading_tr.interval)
    intervals = np.vstack((intervals, [left, right]))
    #nwk_tr_ls = np.append(nwk_tr_ls, leading_tr.newick(precision=5))
    nwk_tr_ls = np.array([leading_tr.newick(precision=5)], dtype=object)

    #SOInvl = left # start of current sliding window; equivalent to intervals[0, 0]
    #EOInvl = right # end of current sliding window; equivalent to intervals[-1, 1]

    for var in ts.variants():
        pos = var.site.position
        gt = var.genotypes
        if np.sum(gt) < min_DAC:
            print(f"{pos}-LF", flush=True)
            continue # ignore low DAF sites
        if pos - intervals[0, 0] < 5e4:
            print(f"{pos}-SKP", flush=True)
            continue # skip beginning of the chromosome

        print(f"{pos}-", flush=True)
        # seek to cover beyond 50,000 bp to the right
        while intervals[-1, 1] - pos < 5e4:
            found = leading_tr.next()
            if found:
                left, right = map(int, leading_tr.interval)
                intervals = np.vstack((intervals, [left, right]))
                nwk_tr_ls = np.append(nwk_tr_ls, leading_tr.newick(precision=5))
                print(f"Seek: {left}-{right}", flush=True)

        if not found: break # reached the end

        # remove trailing trees
        obsolete = (intervals[:, 1] < pos - 5e4)
        if np.sum(obsolete) > 0:
            print(f"Drop: {np.sum(obsolete)} intvls, last-{intervals[obsolete][-1]}", flush=True)
            intervals = intervals[np.logical_not(obsolete)]
            nwk_tr_ls = nwk_tr_ls[np.logical_not(obsolete)]

        fea_at_var = utils.ARG2feature_win(intervals, nwk_tr_ls, pos, gt, win_l, win_r, discreT)
        lin_df = np.concatenate((lin_df, fea_at_var[np.newaxis]))
        pos_ls.append(pos)
        daf_ls.append(np.mean(gt))
        print(gt, lin_df.tolist(), sep='\n', flush=True)

    #np.savez_compressed(out_pref, POS=pos_ls, DAF=daf_ls, FEA=lin_df)
    print("<>", len(daf_ls), len(pos_ls), lin_df.shape)
    
    return 0

sys.exit(main(sys.argv))
