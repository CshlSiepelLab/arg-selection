#!/usr/bin/env python3

import os
import subprocess
import numpy as np
import dendropy
import tskit

RELATE_PATH = '/sonas-hs/siepel/hpc_norepl/home/mo/sel_coef_empirical/relate_v1.0.16_x86_64_static/bin/'
#RELATE_PATH = 'relate_v1.0.16_MacOSX/bin/'
mut_rate = "2.5e-8"
discretT = np.loadtxt('time.txt')
discretT = discretT.astype(int)

def run_RELATE(pp, gtm, Ne): # pp- physical position; gtm: genotype matrix
    # create RELATE input files

    with open("temp.haps", 'w') as hapF:
        for i in range (len(pp)):
            print(str(1), "SNP"+str(i+1), int(pp[i]), "A", "T", *gtm[i], sep=" ", file=hapF)
    # Unpacking Argument Lists: *args

    no_dips = round(gtm.shape[1]/2)
    with open("temp.sample", 'w') as sampF:
        print('ID_1', 'ID_2', 'missing', file=sampF)
        print('0', '0', '0', file=sampF)
        for i in range(no_dips):
            print('UNR'+str(i+1), 'UNR'+str(i+1), 0, file=sampF)

    rho_cMpMb = 1.25
    with open("temp.map", 'w') as mapF:
        ppos = 0
        rdist = 0
        print('pos', 'COMBINED_rate', 'Genetic_map', file=mapF)
        print(int(ppos), rho_cMpMb, rdist, file=mapF)
        
        for loc in pp:
            next_ppos = loc
            rdist = rdist + rho_cMpMb/1e6*(next_ppos - ppos)
            ppos = next_ppos
            print(int(ppos), rho_cMpMb, rdist, file=mapF)

    # run RELATE

    cmd = [RELATE_PATH+"Relate", "--mode", "All",
            "-m", mut_rate,
            "-N", Ne,
            "--haps", "temp.haps",
            "--sample", "temp.sample",
            "--map", "temp.map",
            "-o", "temp"]
    relate_proc = subprocess.run(cmd)
    relate_proc.check_returncode()
    #if relate_proc.returncode != 0:
    #    raise RuntimeError("Relate failed with Exit Stat", relate_proc.returncode)

    # convert to tree-sequence

    cmd = [RELATE_PATH+"RelateFileFormats", "--mode", "ConvertToTreeSequence",
            "-i", "temp",
            "-o", "temp"]
    conv_proc = subprocess.run(cmd)
    conv_proc.check_returncode()

    # load tree sequence

    ts_inferred = tskit.load("temp.trees")

    # clean-up
    os.remove("temp.haps")
    os.remove("temp.sample")
    os.remove("temp.map")
    os.remove("temp.anc")
    os.remove("temp.mut")
    os.remove("temp.trees")

    return ts_inferred

def xtract_fea(tree, var, base):

    # delta= 0.001 
    # tmax = 20000
    # K = 1000
    # discretT = []
    # for i in range(2,K+2):
    #     discretT.append((np.exp(i/(K-1)*np.log(1+delta*tmax))-1)/delta)
    # discretT = np.round(discretT)

    distance_from_root = np.max(tree.calc_node_root_distances()) - discretT
    ones = np.where(var==1)[0] + base
    Canc = []
    Cder = []
    for d in distance_from_root:
        num_lineages_anc = 0
        num_lineages_der = 0
        for nd in tree.preorder_node_iter():
            if not nd._parent_node:
                pass
            else:
                leaf_labels = []
                if nd.distance_from_root() == d:
                    for lnode in nd.leaf_nodes():
                        leaf_labels.append(int(lnode.taxon.label))
                    elem = np.isin(leaf_labels,ones)
                    if len(np.where(elem==True)[0])==len(elem):
                        num_lineages_der += 1
                    else:
                        num_lineages_anc += 1
                elif nd.distance_from_root() >= d and nd._parent_node.distance_from_root() < d:
                    for lnode in nd.leaf_nodes():
                        leaf_labels.append(int(lnode.taxon.label))
                    elem = np.isin(leaf_labels,ones)
                    if len(np.where(elem==True)[0])==len(elem):
                        num_lineages_der += 1
                    else:
                        num_lineages_anc += 1
        Canc.append(num_lineages_anc)
        Cder.append(num_lineages_der)

    C = np.vstack((Canc, Cder))

    return C

def infer_ARG_fea(pos_ls, geno_mtx, put_sel_var, Ne):
    '''Format input, run RELATE on variants of a simulated region and extract features of a region from inferred tree sequence'''

    # convert coordinate and resolve rounding duplicate
    p = pos_ls * 1e5
    p = np.round(p).astype(int)
    while len(p) != len(np.unique(p)):
        for k in range(1,len(p)):
            if p[k] == p[k-1]:
                p[k] = p[k-1] + 1
                p = np.sort(p)

    ts_inferred = run_RELATE(p, geno_mtx, Ne)

    for tree in ts_inferred.trees():
        left, right = map(int, tree.interval)
        if left <= 50000 and right > 50000:
            t = tree.newick(precision=1)
            break
    dp_tree = dendropy.Tree.get(data=t, schema="newick")

    feature_mtx = xtract_fea(dp_tree, put_sel_var, 1) # 2 x 100

    return feature_mtx

def infer_ARG_fea_3(pos_ls, geno_mtx, put_sel_var, Ne):
    '''Format input, run RELATE on variants of a simulated region and extract features of a region from inferred tree sequence'''

    # convert coordinate and resolve rounding duplicate
    p = pos_ls * 1e5
    p = np.round(p).astype(int)
    while len(p) != len(np.unique(p)):
        for k in range(1,len(p)):
            if p[k] == p[k-1]:
                p[k] = p[k-1] + 1
                p = np.sort(p)

    ts_inferred = run_RELATE(p, geno_mtx, Ne)

    center_found = False
    for tree in ts_inferred.trees():
        left, right = map(int, tree.interval)
        
        if center_found:
            rTree = tree.newick(precision=1)
            break
        
        if left <= 50000 and right > 50000:
            cTree = tree.newick(precision=1)
            center_found = True
            a, b = left, right
            continue
            
        lTree = tree.newick(precision=1)

    l_var_pos = p[p<a].max()
    r_var_pos = p[p>b].min()
    l_var_idx = np.where(p==l_var_pos)[0][0]
    c_var_idx = np.argmin(np.absolute(p-50000))
    r_var_idx = np.where(p==r_var_pos)[0][0]

    l_tree = dendropy.Tree.get(data=lTree, schema="newick")
    c_tree = dendropy.Tree.get(data=cTree, schema="newick")
    r_tree = dendropy.Tree.get(data=rTree, schema="newick")

    l_fea = xtract_fea(l_tree, geno_mtx[l_var_idx, :], 1)
    c_fea = xtract_fea(c_tree, geno_mtx[c_var_idx, :], 1)
    r_fea = xtract_fea(r_tree, geno_mtx[r_var_idx, :], 1)

    surr_fea = (l_fea+r_fea)/2
    feature_mtx = np.vstack((c_fea, surr_fea)) # 4 by 100

    return feature_mtx

def infer_ARG_fea_5(pos_ls, geno_mtx, put_sel_var, Ne):
    '''Format input, run RELATE on variants of a simulated region and extract features of a region from inferred tree sequence
        Features include # of lineages at discretized time points & length of non-recomb. segment of surrounding gene trees,
        as well as the # of anc. & der. lineages at discretized time points, length of n.r.s. and der. allelic freq. at the focal site
    '''

    # convert coordinate and resolve rounding duplicate
    p = pos_ls * 1e5
    p = np.round(p).astype(int)
    while len(p) != len(np.unique(p)):
        for k in range(1,len(p)):
            if p[k] == p[k-1]:
                p[k] = p[k-1] + 1
                p = np.sort(p)

    ts_inferred = run_RELATE(p, geno_mtx, Ne)

    trees = []
    intervals = np.empty((0,2), int)
    for tree in ts_inferred.trees():
        left, right = map(int, tree.interval)
        intervals = np.vstack((intervals, [left, right]))
        trees = np.append(trees, tree.newick(precision=1))

    c_idx = np.where((intervals[:,0]<=50000) & (intervals[:,1]>50000))[0][0]
    c_tree = dendropy.Tree.get(data=trees[c_idx], schema="newick")
    c_fea = xtract_fea(c_tree, put_sel_var, 1)

    feature_mtx_5 = np.hstack((c_fea, np.reshape([intervals[c_idx, 1]-intervals[c_idx, 0], np.sum(put_sel_var)/np.shape(put_sel_var)[0]], (2,1))))
    #print(feature_mtx_5.shape)
    s_indices = [c_idx-2, c_idx-1, c_idx+1, c_idx+2]

    for st_idx in s_indices:
        if st_idx < 0 or st_idx >= trees.shape[0]: # in case there are fewer than 5 gene trees in the region
            feature_mtx_5 = np.vstack((feature_mtx_5, np.append(np.sum(c_fea, axis=0), 0)))
            continue
        st = dendropy.Tree.get(data=trees[st_idx], schema="newick")
        root_h = st.max_distance_from_root()
        T = root_h - discretT
        st_fea = np.array([st.num_lineages_at(t) for t in T]+[intervals[st_idx, 1]-intervals[st_idx, 0]])
        feature_mtx_5 = np.vstack((feature_mtx_5, st_fea))
        #print(feature_mtx_5.shape)

    return feature_mtx_5    # 6 x 101