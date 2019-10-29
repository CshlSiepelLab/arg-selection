#!/usr/bin/env python3

import os
import subprocess
import numpy as np
#import pandas as pd
import allel
import dendropy
import tskit

RELATE_PATH = '/sonas-hs/siepel/hpc_norepl/home/mo/sel_coef_empirical/relate_v1.0.16_x86_64_static/bin/'
#RELATE_PATH = '/Users/mo/Google Drive/Cloud Literally/Late_2019/inf_gt/relate_v1.0.16_MacOSX/bin/'
mut_rate = "2.5e-8"

#discretT = np.loadtxt('time.txt')
#discretT = discretT.astype(int)
#K = len(discretT)
delta= 0.001 
tmax = 20000
K = 1000
discretT = []
for i in range(2,K+2):
    discretT.append((np.exp(i/(K-1)*np.log(1+delta*tmax))-1)/delta)
discretT = np.round(discretT)

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

def calc_H1(gt_mtx, pos):

    pos = np.around(pos * 100000) # convert position to coordinate in 100kb region
    hArr = allel.HaplotypeArray(gt_mtx)
    acArr = hArr.count_alleles() 

    H1, H12, H123, H2H1 = allel.garud_h(hArr)

    return H1

def infer_ARG_fea(pos_ls, geno_mtx, put_sel_var, var_pos, Ne, no_ft, norm_iHS):
    '''Format input, run RELATE on variants of a simulated region and extract features of a region from inferred tree sequence
        Features include # of lineages at discretized time points & length of non-recomb. segment of surrounding gene trees,
        as well as the # of anc. & der. lineages at discretized time points, length of n.r.s. and der. allelic freq. at the focal site

        For neutral regions, `var_pos` is the site chosen to have a matching allele frequency;
        for sweep regions, `var_pos` is specified to be 0.5
        `no_ft` - number of flanking gene trees to include on EACH side for feature extraction
        `norm_iHS` - list of normalized iHS score at each variant, cols=["pos", "iHS"]
    '''

    # convert coordinate and resolve rounding duplicate
    p = pos_ls * 1e5
    p = np.round(p).astype(int)
    while len(p) != len(np.unique(p)):
        for k in range(1,len(p)):
            if p[k] == p[k-1]:
                p[k] = p[k-1] + 1
                p = np.sort(p)

    if var_pos == 0.5:
        var_ppos = 50000
    else:
        var_idx = np.where(pos_ls == var_pos)[0][0]
        var_ppos = p[var_idx]

    ts_inferred = run_RELATE(p, geno_mtx, Ne)

    trees = []
    intervals = np.empty((0,2), int)
    for tree in ts_inferred.trees():
        left, right = map(int, tree.interval)
        intervals = np.vstack((intervals, [left, right]))
        trees = np.append(trees, tree.newick(precision=1))

    lin_mtx = np.empty((0, K))
    stat_mtx = np.empty((0, 4))

    c_idx = np.where((intervals[:,0]<=var_ppos) & (intervals[:,1]>var_ppos))[0][0]
    window = range(c_idx-no_ft, c_idx+no_ft+1)
    s_indices = np.take(np.arange(trees.shape[0]), window, mode='clip') # mode='clip' or 'wrap'

    DAF_ls = np.mean(geno_mtx, axis=1)
    pos_iHS = np.round(norm_iHS[:, 0]*1e5)

    for cnt, st_idx in enumerate(s_indices):
        st = dendropy.Tree.get(data=trees[st_idx], schema="newick")
        end = intervals[st_idx, 1]
        begin = intervals[st_idx, 0]

        length = end - begin
        iHS = np.mean(norm_iHS[(pos_iHS>=begin) & (pos_iHS<end), 1])
        H1 = calc_H1(geno_mtx[(p>=begin) & (p<end), :], pos_ls[(p>=begin) & (p<end)])

        if st_idx == c_idx and cnt == no_ft:
            DAF = np.sum(put_sel_var)/np.shape(put_sel_var)[0]
            c_fea = xtract_fea(st, put_sel_var, 1)
            lin_mtx = np.vstack((lin_mtx, c_fea))
            stat_mtx = np.vstack((stat_mtx, [length, iHS, H1, DAF], [length, iHS, H1, DAF]))
        else:
            DAF = np.mean(DAF_ls[(p>=begin) & (p<end)])
            root_h = st.max_distance_from_root()
            T = root_h - discretT
            st_fea = np.array([st.num_lineages_at(t) for t in T])
            lin_mtx = np.vstack((lin_mtx, st_fea))
            stat_mtx = np.vstack((stat_mtx, [length, iHS, H1, DAF]))

    return lin_mtx, stat_mtx    # dim(lin_mtx)=(2*no_ft+2, K); dim(stat_mtx)=(2*no_ft+2, 4)
