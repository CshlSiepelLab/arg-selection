#!/usr/bin/env python3

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "..")) 

import dendropy
import tskit
import pickle
import gzip
import numpy as np

import feature_util as fea

helpMsg = '''
        usage: $./genFeatures.py <.trees PATH> <out prefix> <no_ft>
'''

def main(args):
    if len(args) != 4:    #3 arguments
        return helpMsg

    ts_file=args[1]
    out_pref=args[2]
    no_ft=args[3]
    
    ts = tskit.load(ts_file)
    
    # feature matrices for variant sites
    pos_ls = []
    lin_ls = []
    stat_ls = []

    # features of trees
    #foc_fea_ls = [] # this is differnet for each variant

    # cache
    ppos_ls = np.array([])
    gt_mtx = np.empty((0, len(ts.samples()))) # of row matches len(ppos_ls)

    intervals = np.empty((0, 2), int) # of row = 2*no_ft + 1
    dp_tr_ls = []  # list of dendropy objects len(dp_tr_ls) = 2*no_ft + 1
    flk_fea_ls = np.empty((0, fea.K)) # flk_fea_ls.shape = (2*no_ft + 1, K)
    #lHaD_ls = [] # np.array([len, H1, avgDAF])
    EOInvl = 0 # end of current interval

    for var in ts.variants():
        pos = var.site.position
        gt = var.genotypes

        if pos < EOInvl:
            ppos_ls = np.append(ppos_ls, pos)
            gt_mtx = np.vstack((gt_mtx, gt))
            continue

        if len(dp_tr_ls) == 2*no_ft+1:
            # resolve the middle segment of the previous window
            pos_vs, lin_vs, stat_vs = fea.vars_ARG_fea(ppos_ls, gt_mtx, intervals, dp_tr_ls, flk_fea_ls, no_ft)
            pos_ls += pos_vs
            lin_ls += lin_vs
            stat_ls += stat_vs

            # Remove leftmost tree from cache
            keep = (ppos_ls >= intervals[1, 0])
            ppos_ls = ppos_ls[keep]
            gt_mtx = gt_mtx[keep, :]

            intervals = intervals[1:, :]
            dp_tr_ls = dp_tr_ls[1:]
            flk_fea_ls = flk_fea_ls[1:, :]

        # variant in a new segment
        tree = ts.at(pos)
        left, right = map(int, tree.interval)

        # populate cache
        ppos_ls = np.append(ppos_ls, pos)
        gt_mtx = np.vstack((gt_mtx, gt))

        intervals = np.vstack((intervals, [left, right]))

        st = dendropy.Tree.get(data=tree, schema="newick")
        root_h = st.max_distance_from_root()
        T = root_h - fea.discretT
        st_fea = np.array([st.num_lineages_at(t) for t in T])

        dp_tr_ls.append(st)
        flk_fea_ls = np.vstack(flk_fea_ls, st_fea)

        EOInvl = right

    with open(out_pref+'.pickle', 'wb') as f:
        pickle.dump((pos_ls, lin_ls, stat_ls), f, pickle.HIGHEST_PROTOCOL)
    
    return 0

sys.exit(main(sys.argv))
