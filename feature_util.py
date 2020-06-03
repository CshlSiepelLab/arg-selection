#!/usr/bin/env python3

import os
import subprocess
import numpy as np
#import pandas as pd
import allel
import dendropy
import tskit

RELATE_PATH = '/sonas-hs/siepel/hpc_norepl/home/mo/relate_v1.0.17_x86_64_static/'
#RELATE_PATH = '~/relate_v1.0.16_MacOSX/'
ii32MAX = np.iinfo(np.int32).max

time_file_path = os.path.dirname(os.path.abspath(__file__))+'/sim2args/time.txt'
discretT = np.loadtxt(time_file_path)
discretT = discretT.astype(int)
K = len(discretT)

# delta= 0.001
# tmax = 20000
# K = 100
# discretT = []
# for i in range(2,K+2):
#     discretT.append((np.exp(i/(K-1)*np.log(1+delta*tmax))-1)/delta)
# discretT = np.round(discretT)

def run_RELATE(pp, gtm, Ne, var_ppos=-1, rho_cMpMb=1.25, mut_rate="2.5e-8"):
    '''pp- physical position; gtm: genotype matrix; var_ppos: physical position of the variant
        var_ppos = -1 : default to skip RELATE selection inference, returns pval=0.1
        rho_cMpMb, mut_rate: default to human
    '''
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

    with open("temp.poplabels", 'w') as popF:
        print('sample', 'population', 'group', 'sex', file=popF)
        for i in range(no_dips):
            print('UNR'+str(i+1), 'CEU', 'EUR', 1, file=popF)

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

    relate_cmd = [RELATE_PATH+"bin/Relate", "--mode", "All",
            "-m", mut_rate,
            "-N", Ne,
            "--haps", "temp.haps",
            "--sample", "temp.sample",
            "--map", "temp.map",
            "-o", "temp",
            #"--memory", "8", # default is 5 GB
            "--seed", "1"]

    popsize_cmd = [RELATE_PATH+"scripts/EstimatePopulationSize/EstimatePopulationSize.sh",
            "-i", "temp",
            "-m", mut_rate,
            "--poplabels", "temp.poplabels",
            "--threshold", "10",
            #"--num_iter", "10", # default is 5
            "-o", "temp_popsize",
            "--seed", "1"]

    wg_cmd = [RELATE_PATH+"scripts/EstimatePopulationSize/EstimatePopulationSize.sh",
            "-i", "temp",
            "-m", mut_rate,
            "--poplabels", "temp.poplabels",
            "--threshold", "0",
            "--coal", "temp_popsize.coal",
            "--num_iter", "1",
            "-o", "temp_wg"]

    conv_cmd = [RELATE_PATH+"bin/RelateFileFormats", "--mode", "ConvertToTreeSequence",
            "-i", "temp_wg",
            "-o", "temp_wg"]

    #clean_cmd = [RELATE_PATH+"bin/Relate", "--mode", "Clean", "-o", "temp"]

    loop_cnt = 0
    # run RELATE
    while True:
        loop_cnt += 1
        if loop_cnt > 20:
            print("Milestone: attempt exceeds 20, ABORTS!!")
            subprocess.call("rm temp*", shell=True)
            return None, None

        print("Milestone: Running RELATE pipeline, try_", loop_cnt, sep='', flush=True)
        relate_cmd[-1] = str(np.random.randint(ii32MAX))
        popsize_cmd[-1] = str(np.random.randint(ii32MAX))
        relate_proc = subprocess.run(relate_cmd)
        if relate_proc.returncode != 0: continue
        #relate_proc.check_returncode()
        # while relate_proc.returncode != 0:
        # #    raise RuntimeError("Relate failed with Exit Stat", relate_proc.returncode)
        #     print("relate_proc_failed, restart")
        #     relate_proc = subprocess.run(cmd)

        # re-estimate branch lengths
        popsize_proc = subprocess.run(popsize_cmd)
        if popsize_proc.returncode != 0: continue
        # output: _popsize.pdf, _popsize.anc.gz, _popsize.mut.gz, _popsize.dist, _popsize.coal. _popsize.bin, _popsize_avg.rate

        # re-estimate branch length for ENTIRE genealogy
        wg_proc = subprocess.run(wg_cmd)
        if wg_proc.returncode != 0: continue

        # convert to tree-sequence
        conv_proc = subprocess.run(conv_cmd)
        if conv_proc.returncode != 0: continue
        # if the conversion throws <time[parent] must be greater than time[child]> error, rerun the last resampling step
        print("Milestone: RELATE pipeline success", flush=True)
        break

    # run RELATE selection inference
    if var_ppos in pp:
        pval = RELATE_sel_inf(var_ppos, mut_rate)
    else:
        pval = 0.1 # selected variant not in the list

    # load tree sequence
    ts_inferred = tskit.load("temp_wg.trees")

    # clean-up
    subprocess.call("rm temp*", shell=True)

    return ts_inferred, pval

def RELATE_sel_inf(locOI, mut_rate):

    # cmd = [RELATE_PATH+"bin/RelateCoalescentRate",
    #         "--mode", "EstimatePopulationSize",
    #         "-i", "temp",
    #         "-o", "temp"]

    # pwcoal_proc = subprocess.run(cmd)
    # pwcoal_proc.check_returncode()
    # # output: temp.coal, temp.bin

    cmd = [RELATE_PATH+"scripts/DetectSelection/DetectSelection.sh",
            "-i", "temp_wg",
            "-o", "temp_sel",
            "-m", mut_rate,
            "--first_bp", str(locOI),
            "--last_bp", str(locOI)]
            #"--coal", "temp.coal"]

    sele_proc = subprocess.run(cmd)
    #sele_proc.check_returncode()
    while sele_proc.returncode != 0:
        print("sele_proc_failed, restart")
        sele_proc = subprocess.run(cmd)
    # output: temp_sel.lin, temp_sel.freq, temp_sel.sele, temp_sel.anc, temp_sel.mut

    sel_pval = np.genfromtxt("temp_sel.sele", skip_header=1)

    if sel_pval.size == 0:
        return 0 # RELATE cannot test for selection at site where the mutation cannot be mapped to a unique branch
    else:
        return sel_pval[-1]

def xtract_fea(tree, var, base):
    '''
    discoal tree tips range from 0 to 197 (base = 0)
    RELATE tree tips range from 1 to 198 (base = 1)
    '''
    # lf_dist = np.floor(tree.calc_node_root_distances()) # round down!
    # if (np.max(lf_dist) - np.min(lf_dist)) > 1:
    #     print("Anomaly: leaf age -", lf_dist, flush=True)
    # distance_from_root = np.min(lf_dist) - discretT
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

def calc_H1(gt_mtx):

    #pos = np.around(pos * 100000) # convert position to coordinate in 100kb region
    hArr = allel.HaplotypeArray(gt_mtx)
    acArr = hArr.count_alleles() 

    H1, H12, H123, H2H1 = allel.garud_h(hArr)

    return H1

def infer_ARG_fea(pos_ls, geno_mtx, put_sel_var, var_pos, Ne, mu_string, rho_1e8, no_ft, norm_iHS):
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

    ts_inferred, RELATE_pval = run_RELATE(p, geno_mtx, Ne, var_ppos, rho_1e8, mu_string)
    if ts_inferred is None:
        return None, None, None

    trees = []
    intervals = np.empty((0,2), int)
    for tree in ts_inferred.trees():
        left, right = map(int, tree.interval)
        intervals = np.vstack((intervals, [left, right]))
        trees = np.append(trees, tree.newick(precision=1))

    intervals[-1, 1] = 1e5 # force last tree to cover the rest of the region

    lin_mtx = np.empty((0, K))
    #stat_mtx = np.empty((0, 4))
    stat_mtx = np.empty((0, 3))
    # iHS is not included for now. iHS is not calculated for sites with minor allele freq. <.05, which creates NaNs in the feature matrix

    c_idx = np.where((intervals[:,0]<=var_ppos) & (intervals[:,1]>var_ppos))[0][0]
    window = range(c_idx-no_ft, c_idx+no_ft+1)
    s_indices = np.take(np.arange(trees.shape[0]), window, mode='clip') # mode='clip' or 'wrap'

    DAF_ls = np.mean(geno_mtx, axis=1)
    #pos_iHS = np.round(norm_iHS[:, 0]*1e5)

    for cnt, st_idx in enumerate(s_indices):
        st = dendropy.Tree.get(data=trees[st_idx], schema="newick")
        end = intervals[st_idx, 1]
        begin = intervals[st_idx, 0]

        length = end - begin
        #iHS = np.mean(norm_iHS[(pos_iHS>=begin) & (pos_iHS<end), 1])
        H1 = calc_H1(geno_mtx[(p>=begin) & (p<end), :])
        avgDAF = np.mean(DAF_ls[(p>=begin) & (p<end)])

        if st_idx == c_idx and cnt == no_ft:
            DAF = np.sum(put_sel_var)/np.shape(put_sel_var)[0]
            c_fea = xtract_fea(st, put_sel_var, 1)
            lin_mtx = np.vstack((lin_mtx, c_fea)) # C = np.vstack((Canc, Cder))
            #stat_mtx = np.vstack((stat_mtx, [length, iHS, H1, avgDAF], [length, iHS, H1, DAF]))
            stat_mtx = np.vstack((stat_mtx, [length, H1, avgDAF], [length, H1, DAF]))
        else:
            # lfrt_dist = np.floor(st.calc_node_root_distances()) # round down!
            # if (np.max(lfrt_dist) - np.min(lfrt_dist)) > 1:
            #     print("Anomaly: leaf age -", lfrt_dist, flush=True)
            # root_h = np.min(lfrt_dist)
            root_h = st.max_distance_from_root()
            T = root_h - discretT
            st_fea = np.array([st.num_lineages_at(t) for t in T])
            lin_mtx = np.vstack((lin_mtx, st_fea))
            #stat_mtx = np.vstack((stat_mtx, [length, iHS, H1, avgDAF]))
            stat_mtx = np.vstack((stat_mtx, [length, H1, avgDAF]))

    return lin_mtx, stat_mtx, RELATE_pval    # dim(lin_mtx)=(2*no_ft+2, K); dim(stat_mtx)=(2*no_ft+2, 4)

def vars_ARG_fea(ppos_ls, gtm, intvls, dp_tr_ls, flk_fea_ls, no_ft, minDAC):
    '''
    ppos_ls: list of physical positions of the variants
    gtm: genotype matrix of the variants
    intvls: intervals of the trees
    dp_tr_ls: dendropy objects of the trees
    flk_fea_ls: precomputed features of the trees (if they are the flanking trees)
    no_ft: # of flanking trees on EACH side
    min_DAC: minimum derived allele count for feature calculation
    '''

    pos_vars = []
    lin_vars = []
    stat_vars = []
    tree_vars = []

    stat_cache = np.empty((0, 3))

    for idx in range(len(dp_tr_ls)):
        end = intvls[idx, 1]
        begin = intvls[idx, 0]
        length = end - begin
        gtm_loc = gtm[(ppos_ls>=begin) & (ppos_ls<end), :]
        H1 = calc_H1(gtm_loc)
        avgDAF = np.mean(gtm_loc)
        stat_cache = np.vstack((stat_cache, [length, H1, avgDAF]))

    foc_sites = (ppos_ls>=intvls[no_ft, 0]) & (ppos_ls<intvls[no_ft, 1])
    focT_nwk = dp_tr_ls[no_ft].as_string(schema="newick")

    for var_i in np.argwhere(foc_sites).flatten(): # argwhere returns a column vector, needs to be flattened!
        gt = gtm[var_i, :]
        if np.sum(gt) < minDAC: continue
        
        DAF = np.mean(gt)
        c_fea = xtract_fea(dp_tr_ls[no_ft], gt, 1)

        lin_mtx = np.vstack((flk_fea_ls[:no_ft, :], c_fea, flk_fea_ls[-no_ft:, :])) # C = np.vstack((Canc, Cder))
        stat_mtx = np.vstack((stat_cache[:no_ft+1, :], stat_cache[no_ft:, :]))
        stat_mtx[no_ft+1, 2] = DAF

        pos_vars.append(ppos_ls[var_i])
        lin_vars.append(lin_mtx)
        stat_vars.append(stat_mtx)
        tree_vars.append(focT_nwk)

    return pos_vars, lin_vars, stat_vars, tree_vars
