#!/usr/bin/env python3

import os
import sys
import dendropy
import numpy as np

helpMsg = '''
    usage: $./discoalT2feature.py <POP> <out_Pref>
'''

def parse_tree_line(line, t_max, delta=1e3, K=50):
    
    intv, nwk = line.split(']')
    intv = int(intv[1:])
    tr = dendropy.Tree.get(data=nwk, schema="newick")

    I = np.arange(K)
    T = (np.exp(I/(K-1)*np.log(1+delta*t_max))-1)/delta
    
    root_h = tr.max_distance_from_root()
    T = root_h - T
    lin = np.array([tr.num_lineages_at(t) for t in T])
    lin[lin==0] = 1 # num_lineages_at(t) outputs zero when t<=0, fix this to 1
    
    return intv, lin, root_h

def parse_hdswp_cmd(cmd_line):
    
    args = cmd_line.split()
    N_curr = int(args[args.index('-N') + 1])
    alpha = float(args[args.index('-a') + 1])
    #AF = float(args[args.index('-c') + 1])
    selcoef = alpha/(2*N_curr)
    
    return selcoef

def process_discoalT(discoal_file, cat, region_size=1e5, win_size=2e4, K=50, t_max=0.5):
    with open(discoal_file, "r") as discoalF:
        if cat == 'swp': 
            sc = parse_hdswp_cmd(discoalF.readline().strip())
        elif cat == 'neu':
            sc = 0

        win_ends = np.arange(win_size-1, region_size, win_size) # 0-based indexing
        region_fea = np.empty((len(win_ends), K)) # feature of the entire region

        ptr = 0 # 0-based indexing
        win=0
        seg_fea = np.empty((0, K)) # (K time points) temporary storage of features of the trees in one segment/window
        seg_size = []

        for ln in discoalF:
            if ln[0] != '[': continue
            intv, fea, ht = parse_tree_line(ln.strip(), t_max=t_max, delta=1e3, K=K)
            left = ptr
            right = ptr + intv
            ptr += intv
            
            if right > win_ends[win]:
                seg_size.append(win_ends[win]-left+1)
                seg_fea = np.vstack((seg_fea, fea))
                region_fea[win, :] = np.average(seg_fea, axis=0, weights=seg_size)
                seg_fea = np.empty((0, K))
                seg_fea = np.vstack((seg_fea, fea))
                seg_size = []
                seg_size.append(right-win_ends[win]+1)
                win += 1
                continue
                
            seg_size.append(right-left)
            seg_fea = np.vstack((seg_fea, fea))
    
    return sc, region_fea

def main(args):
    if len(args) != 3:    #2 argument
        return helpMsg
    
    pop_pref=args[1] # e.g. CHB/CHB18596
    
    dim = 5*50 # no. of windows x no. of time points
    
    swp_samp = 5000
    X_swp = np.zeros((swp_samp, dim))
    y_swp = np.empty(swp_samp)
    
    for samp in range(swp_samp):
        discoalF_path = pop_pref+'_hdswp_'+str(samp+1)+'.discoal'
        if os.stat(discoalF_path).st_size == 0: continue
        sc, fea = process_discoalT(discoalF_path, 'swp')
        
        X_swp[samp, :] = fea.flatten()
        y_swp[samp] = sc
    
    neu_samp = 5000
    X_neu = np.empty((neu_samp, dim))
    y_neu = np.zeros(neu_samp)

    for samp in range(neu_samp):
        discoalF_path = pop_pref+'_neutral_'+str(samp+1)+'.discoal'
        sc, fea = process_discoalT(discoalF_path, 'neu')
        
        X_neu[samp, :] = fea.flatten()
        y_neu[samp] = sc
        
    outPref=args[2]
    np.savez_compressed(outPref, X_swp=X_swp, y_swp=y_swp, X_neu=X_neu, y_neu=y_neu)

sys.exit(main(sys.argv))
