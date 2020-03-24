#!/usr/bin/env python3

import numpy as np
#import pandas as pd
import dendropy
import tskit

def get_site_ppos(ts):
    var_ppos_ls = []
    prev_pos = 0
    for site in ts.sites():
        site_pos = int(site.position)
        if site_pos <= prev_pos:
            if prev_pos == 49999:
                var_ppos_ls.append(-1) # flag indicating this site should be removed
                continue
            else:
                site_pos = prev_pos + 1
        var_ppos_ls.append(site_pos)
        prev_pos = site_pos
    return np.array(var_ppos_ls)

def samp_var(geno_mtx, pos_ls, AF_min=0.2, AF_max=0.9):
    AF = np.mean(geno_mtx, axis=1)
    target = np.random.uniform(AF_min, AF_max)
    for idx in np.argsort(np.abs(AF - target)):
        if pos_ls[idx] > 45000 and pos_ls[idx] < 55000:
            return idx

def discretizeT(delta=0.01, tmax=1e5, K=100):
    discretT = []
    for i in range(K):
        discretT.append((np.exp(i/(K-1)*np.log(1+delta*tmax))-1)/delta)
    return np.round(discretT)

def cnt_anc_der_lin(tree, var, time_pts):
    lf_dist = np.floor(tree.calc_node_root_distances()) # round down!
    if (np.max(lf_dist) - np.min(lf_dist)) > 1:
        print("Anomaly: leaf age -", lf_dist, flush=True)
    distance_from_root = np.min(lf_dist) - time_pts
    taxOI = np.where(var==1)[0] + 1 # By default, leaf nodes are labelled with their numerical ID + 1
    Canc = []
    Cder = []
    Cmix = []
    for d in distance_from_root:
        num_lineages_anc = 0
        num_lineages_der = 0
        num_lineages_mix = 0
        for nd in tree.preorder_node_iter():
            if not nd._parent_node: continue # root node
            edge_cross_time = (nd.root_distance == d or (nd.root_distance > d and nd._parent_node.root_distance < d))
            if edge_cross_time:
                leaves = []
                for lnode in nd.leaf_iter(): leaves.append(int(lnode.taxon.label))
                der_cnt = np.isin(leaves, taxOI)
                if np.sum(der_cnt) == der_cnt.shape[0]:
                    num_lineages_der += 1
                elif np.sum(der_cnt) == 0:
                    num_lineages_anc += 1
                else:
                    num_lineages_mix += 1
                    
        Canc.append(num_lineages_anc)
        Cder.append(num_lineages_der)
        Cmix.append(num_lineages_mix)
        
    return np.vstack((Canc, Cmix, Cder))

def ts2feature(ts, vOI_pos, vOI_gt, ws_offset, we_offset, time_pts):
    K = len(time_pts)
    W = len(ws_offset)
    region_fea = np.empty((W*3, K)) # (Canc, Cmix, Cder) for each window
    
    win_start = (vOI_pos + ws_offset).astype(int)
    win_end = (vOI_pos + we_offset).astype(int)
    win = 0
    seg_fea = np.empty((0, 3, K)) # temp storage of features of the trees in one segment/window
    seg_size = []
    
    for tree in ts.trees():
        left, right = map(int, tree.interval) # apply function int() to tree.interval, left is inclusive and right is exclusive
        if right <= win_start[0]: continue
        dp_tree = dendropy.Tree.get(data=tree.newick(precision=2), schema="newick")
        fea = cnt_anc_der_lin(dp_tree, vOI_gt, time_pts) # shape: (3, K)
        
        while right >= win_end[win]:
            seg_size.append(win_end[win]-max(left, win_start[win]))
            seg_fea = np.concatenate((seg_fea, np.array([fea])))
            
            if seg_fea.shape[0] > 1:
                region_fea[win*3:win*3+3, :] = np.average(seg_fea, axis=0, weights=seg_size)
            else:
                region_fea[win*3:win*3+3, :] = seg_fea # broadcasting
            
            win += 1
            if win == W: return region_fea
            seg_fea = np.empty((0, 3, K))
            seg_size = []
        
        intvl = right-max(left, win_start[win])
        if intvl == 0: continue
        seg_size.append(intvl)
        seg_fea = np.concatenate((seg_fea, np.array([fea])))