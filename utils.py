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

def discrete_pos(pos_ls, reg_len=1e5):
    ''' convert coordinate and resolve rounding duplicate
    pos_ls should be (0, 1)
    TO BE TESTED
    '''
    ppos = pos_ls*reg_len
    ppos = np.round(ppos).astype(int)

    for p in range(1, len(ppos)):
        if ppos[p] <= ppos[p-1]:
            ppos[p] = ppos[p-1] + 1

    return ppos

def samp_var(geno_mtx, pos_ls, AF_min=0.2, AF_max=0.9, left=46000, right=54000):
    AF = np.mean(geno_mtx, axis=1)
    target = np.random.uniform(AF_min, AF_max)
    for idx in np.argsort(np.abs(AF - target)):
        if pos_ls[idx] > left and pos_ls[idx] < right: # make sure the ends are covered by gene trees
            return idx

def discretizeT(delta=0.01, tmax=1e5, K=100):
    discretT = []
    for i in range(K):
        discretT.append((np.exp(i/(K-1)*np.log(1+delta*tmax))-1)/delta)
    return np.round(discretT)

def cnt_anc_der_lin(tree, var, time_pts, mix=True, base=1):
    '''
    ## LATEST feature extraction method ##
    discoal tree tips range from 0 to 197 (base = 0)
    RELATE tree tips range from 1 to 198 (base = 1)
    '''
    lf_dist = np.floor(tree.calc_node_root_distances()) # round down!
    if (np.max(lf_dist) - np.min(lf_dist)) > 1:
        print("Anomaly: leaf age -", lf_dist, flush=True)
    distance_from_root = np.min(lf_dist) - time_pts
    taxOI = np.where(var==1)[0] + base # By default, leaf nodes are labelled with their numerical ID + 1
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
        
    if mix:
        return np.vstack((Canc, Cmix, Cder))
    else:
        return np.vstack((Canc, Cder))

def ARG2feature(intervals, trees, var_ppos, vOI_gt, no_ft, time_pts, taxa):
    ''' Flanking tree approach
    intervals: dataframe with start and end position of genealogies
    trees: list of newick strings
    var_ppos: integer physical position of focal variant
    vOI_gt: genotype of focal variant
    no_ft: # of flanking trees for feature extraction
    time_pts: array of discretized time points
    taxa: starting number of taxa (0 for discoal or 1 for RELATE)
    '''

    lin_mtx = np.empty((0, K))
    #stat_mtx = np.empty((0, 3)) # omitted sum stats for now

    c_idx = np.where((intervals[:,0]<=var_ppos) & (intervals[:,1]>var_ppos))[0][0]
    window = range(c_idx-no_ft, c_idx+no_ft+1)
    s_indices = np.take(np.arange(trees.shape[0]), window, mode='clip') # mode='clip' or 'wrap'

    #DAF_ls = np.mean(geno_mtx, axis=1)

    for cnt, st_idx in enumerate(s_indices):
        st = dendropy.Tree.get(data=trees[st_idx], schema="newick")
        end = intervals[st_idx, 1]
        begin = intervals[st_idx, 0]

        # length = end - begin
        # H1 = calc_H1(geno_mtx[(p>=begin) & (p<end), :])
        # avgDAF = np.mean(DAF_ls[(p>=begin) & (p<end)])

        if st_idx == c_idx and cnt == no_ft:
            # DAF = np.sum(put_sel_var)/np.shape(put_sel_var)[0]
            c_fea = cnt_anc_der_lin(st, vOI_gt, time_pts, mix=False, base=taxa)
            lin_mtx = np.vstack((lin_mtx, c_fea)) # C = np.vstack((Canc, Cder))
            #stat_mtx = np.vstack((stat_mtx, [length, H1, avgDAF], [length, H1, DAF]))
        else:
            lfrt_dist = np.floor(st.calc_node_root_distances()) # round down!
            if (np.max(lfrt_dist) - np.min(lfrt_dist)) > 1:
                print("Anomaly: leaf age -", lfrt_dist, flush=True)
            root_h = np.min(lfrt_dist)
            #root_h = st.max_distance_from_root()
            T = root_h - time_pts
            st_fea = np.array([st.num_lineages_at(t) for t in T])
            lin_mtx = np.vstack((lin_mtx, st_fea))
            #stat_mtx = np.vstack((stat_mtx, [length, H1, avgDAF]))

    return lin_mtx # dim(lin_mtx)=(2*no_ft+2, K)

def ts2feature(ts, vOI_pos, vOI_gt, ws_offset, we_offset, time_pts):
    ''' Full region window approach '''
    K = len(time_pts)
    W = len(ws_offset)
    region_fea = np.empty((W*3, K)) # (Canc, Cmix, Cder) for each window
    
    win_start = (vOI_pos + ws_offset).astype(int)
    win_end = (vOI_pos + we_offset).astype(int)
    win = 0
    seg_fea = np.empty((0, 3, K)) # temp storage of features of the trees in one segment/window
    seg_size = []
    
    tcnt = 0
    ttot = ts.num_trees
    for tree in ts.trees():
        tcnt += 1
        left, right = map(int, tree.interval) # apply function int() to tree.interval, left is inclusive and right is exclusive
        if right <= win_start[0]: continue
        if tcnt == ttot: right = 1e5 # coerce the last tree to cover the entire region
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