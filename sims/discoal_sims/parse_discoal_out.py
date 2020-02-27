#!/usr/bin/env python3

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

#import dendropy
import numpy as np
import h5py

import feature_util as fea

helpMsg = '''
    usage: $./parse_discoal_out.py <`-s` or `-n`> <handle> <no_repl> <out_dir> <no_partition> 

    Parse all discoal_<handle>/discoal_<handle>_*.discoal files
    Outputs <out_dir>/discoal_<handle>_tru.hdf5 or .npz (feature from true gene trees)

                SC, CAF, *AFtraj, fea_Mtx [<= NRseg_len, trees]

            <out_dir>/discoal_<handle>_inf.hdf5 or .npz (genotypes for ARG inference)

                SC, CAF, *AFtraj, var_Pos, gt_Mtx, focVar_idx

            * for sweep simulations only

    <no_repl>: Total number of simulation runs
    <np_partition>: Number of partitions
    -s for sweep simulations
    -n for neutral simulations
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

def process_discoalT(discoal_file, cat, K, delta, t_max, no_fl_win, win_size, region_size):

    # Pass 1: get sc, afTraj, var_pos and gt_mtx
    with open(discoal_file, "r") as discoalF:
        if cat == 'swp': 
            sc = parse_hdswp_cmd(discoalF.readline().strip())
        elif cat == 'neu':
            sc = 0

    # Pass 2: Feature extraction
    with open(discoal_file, "r") as discoalF:
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
    
    return SC, CAF, AFtraj, feaMtx, varPos, GtMtx, focVarIdx # AFtraj should be None if in `-n` mode

def main(args):
    if len(args) != 6:    #5 argument
        return helpMsg
    
    mode = args[1]
    handle = args[2]
    no_sims = int(args[3])
    out_dir = args[4]
    no_part = int(args[5])
    
    ## Hyper-parameters ##
    K = 100 # number of time points for feature extraction
    delta = 0.001
    t_max = 0.0005 # in 4N_0 unit
    no_fl_win = 2 # number of flanking windows on each side for fea. extr.
    win_size = 1e3 # size of each window (in bp)
    region_size = 1e5

    SC_arr = np.empty(no_sims) # of floats
    CAF_arr = np.empty(no_sims) # of floats
    #AFtraj can have uneven dimension
    
    ## for tru ##
    feaMtx_arr = np.empty((no_sims, 2*(2*no_fl_win+1), K)) # 3 dimensional

    ## for inf ##
    varPos_arr = np.empty(no_sims) # of floats
    #GtMtx can have uneven dimension
    focVarIdx_arr = np.empty(no_sims) # of ints
    
    #The File object does double duty as the HDF5 root group, and serves as your entry point into the file:
    f_tru = h5py.File(out_dir+'/discoal_'+handle+'_tru.hdf5', "w")
    f_inf = h5py.File(out_dir+'/discoal_'+handle+'_inf.hdf5', "w")

    for samp in range(no_sims):
        discoalF_path = 'discoal_'+handle+'/discoal_'+handle+'_'+str(samp)+'.discoal'
        #if os.stat(discoalF_path).st_size == 0: continue
        SC, CAF, AFtraj, feaMtx, varPos, GtMtx, focVarIdx = process_discoalT(discoalF_path, mode, K, delta, t_max, no_fl_win, win_size, region_size)
        # AFtraj should be None if in `-n` mode

        SC_arr[samp] = SC
        CAF_arr[samp] = CAF
        if mode == '-s':
            f_tru.create_dataset('AF_'+str(samp), data=AFtraj)
            f_inf.create_dataset('AF_'+str(samp), data=AFtraj)

        feaMtx_arr[samp, :, :] = feaMtx

        varPos_arr[samp] = varPos
        #f_tru.create_dataset('gtm_'+str(samp), data=GtMtx)
        f_inf.create_dataset('gtm_'+str(samp), data=GtMtx)
        focVarIdx_arr[samp] = focVarIdx        

        #np.savez_compressed(outPref, X_swp=X_swp, y_swp=y_swp, X_neu=X_neu, y_neu=y_neu)

    #initialize the dataset to an existing NumPy array by providing the data parameter
    f_tru.create_dataset('SC_arr', data=SC_arr)
    f_tru.create_dataset('CAF_arr', data=CAF_arr)
    f_tru.create_dataset('feaMtx_arr', data=feaMtx_arr)

    f_inf.create_dataset('SC_arr', data=SC_arr)
    f_inf.create_dataset('CAF_arr', data=CAF_arr)
    f_inf.create_dataset('varPos_arr', data=varPos_arr)
    f_inf.create_dataset('focVarIdx_arr', data=focVarIdx_arr)

    f_tru.close()
    f_inf.close()

sys.exit(main(sys.argv))
