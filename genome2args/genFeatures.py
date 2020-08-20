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
        usage: $./genFeatures.py <.trees PATH> <.meta PATH> <out prefix> <no_ft> <min_DAF>
'''

time_file_path =  os.path.join(os.path.dirname(__file__), "..", "time.txt")
discretT = np.loadtxt(time_file_path)
discretT = discretT.astype(int)
K = len(discretT)
#no_ft = 2 # of flanking gene trees to include on EACH side for feature extraction

def main(args):
    if len(args) != 6:    #5 arguments
        return helpMsg

    ts_file=args[1]
    meta_file=args[2]
    out_pref=args[3]
    no_ft=int(args[4])
    min_DAF=float(args[5])
    
    ts = tskit.load(ts_file)
    # meta = np.genfromtxt(meta_file, dtype=None, encoding=None) 
    # rs_ID | pos | code e.g. ('rs59993883', 135500045, 0)
    meta_ID = np.loadtxt(meta_file, dtype=str, usecols=0, encoding=None)
    meta_pos = np.loadtxt(meta_file, dtype=int, usecols=1)
    meta_code = np.loadtxt(meta_file, dtype=int, usecols=2)
    print(f"<> No. sites in: ts_file-{ts.num_sites}, meta_file-{meta_ID.shape[0]}; {meta_pos.shape[0]}; {meta_code.shape[0]}", flush=True)

    no_haptypes = ts.num_samples
    min_DAC = round(min_DAF*no_haptypes)
    print(f"<> Min DAC threshold: {min_DAC}/{no_haptypes}")
    
    # feature matrices for variant sites
    id_ls = []
    daf_ls = []
    gt_ls = []
    pos_ls = []
    lin_ls = []
    #stat_ls = []
    tree_ls = [] #nwk strings of the trees at each site

    # features of trees
    #foc_fea_ls = [] # this is differnet for each variant

    # cache
    ppos_ls = np.array([])
    id_cache = np.array([])
    #daf_cache = np.array([])
    code_cache = np.array([])
    gt_mtx = np.empty((0, no_haptypes), int) # of row matches len(ppos_ls)

    intervals = np.empty((0, 2), int) # of row = 2*no_ft + 1
    dp_tr_ls = []  # list of dendropy objects len(dp_tr_ls) = 2*no_ft + 1
    flk_fea_ls = np.empty((0, K)) # flk_fea_ls.shape = (2*no_ft + 1, K)
    #lHaD_ls = [] # np.array([len, H1, avgDAF])
    EOInvl = 0 # end of current interval

    for var in ts.variants():
        pos = var.site.position
        gt = var.genotypes
        if np.sum(gt) == 0: continue # ignore invariant sites

        if pos < EOInvl:
            ppos_ls = np.append(ppos_ls, pos)
            meta_idx = np.nonzero(meta_pos == pos)[0][0]
            id_cache = np.append(id_cache, meta_ID[meta_idx])
            code_cache = np.append(code_cache, meta_code[meta_idx])
            gt_mtx = np.vstack((gt_mtx, gt))
            continue

        if len(dp_tr_ls) == 2*no_ft+1:
            # resolve the middle segment of the previous window
            pos_vs, lin_vs, _, tree_vs, id_vs, daf_vs, gt_vs = utils.vars_ARG_fea(ppos_ls, id_cache, code_cache, gt_mtx, intervals, dp_tr_ls, flk_fea_ls, no_ft, min_DAC, discretT)
            id_ls += id_vs
            daf_ls += daf_vs
            gt_ls += gt_vs
            pos_ls += pos_vs
            lin_ls += lin_vs
            #stat_ls += stat_vs
            tree_ls += tree_vs

            # Remove leftmost tree from cache
            keep = (ppos_ls >= intervals[1, 0])
            ppos_ls = ppos_ls[keep]
            id_cache = id_cache[keep]
            code_cache = code_cache[keep]
            gt_mtx = gt_mtx[keep, :]

            intervals = intervals[1:, :]
            dp_tr_ls = dp_tr_ls[1:]
            flk_fea_ls = flk_fea_ls[1:, :]

        # variant in a new segment
        tree = ts.at(pos)
        left, right = map(int, tree.interval)

        # populate cache
        ppos_ls = np.append(ppos_ls, pos)
        meta_idx = np.nonzero(meta_pos == pos)[0][0]
        id_cache = np.append(id_cache, meta_ID[meta_idx])
        code_cache = np.append(code_cache, meta_code[meta_idx])
        gt_mtx = np.vstack((gt_mtx, gt))

        intervals = np.vstack((intervals, [left, right]))

        st = dendropy.Tree.get(data=tree.newick(precision=1), schema="newick")
        lfrt_dist = st.calc_node_root_distances()
        lfrt_dist = np.floor(lfrt_dist) # round down!
        root_h = np.min(lfrt_dist)
        T = root_h - discretT
        st_fea = np.array([st.num_lineages_at(t) for t in T])

        dp_tr_ls.append(st)
        flk_fea_ls = np.vstack((flk_fea_ls, st_fea))

        EOInvl = right

    with open(out_pref+'.pickle', 'wb') as f:
        pickle.dump((id_ls, daf_ls, pos_ls, lin_ls, tree_ls, gt_ls), f, pickle.HIGHEST_PROTOCOL)

    print("<>", len(id_ls), len(daf_ls), len(pos_ls), len(lin_ls), len(tree_ls), len(gt_ls))
    
    return 0

sys.exit(main(sys.argv))
