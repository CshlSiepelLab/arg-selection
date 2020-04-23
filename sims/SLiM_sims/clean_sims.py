#!/usr/bin/env python3

import numpy as np
import sys, os

helpMsg = '''
        usage: $./clean_sims.py [mode] [handle] [sim_p_thr] [tot_thr] > [logfile.txt]
            Mode:   '[thr]' - clean .trees files of a PARTICULAR thread without a EXITSTAT_0 sanity check
                    '[meta_hdl]' - resolve invalid `%%` lines in [meta_hdl].rawmeta
                                   outputs [meta_hdl].meta
            Make sure [handle].0exit contains recap_EXITSTAT_0 for ALL successful sims
'''

def clean_meta_lines(met_hdl, sarr):
    meta_kept = 0
    meta_removed = 0
    with open(met_hdl+".rawmeta", 'r') as rmF:
        with open(met_hdl+".meta", 'w') as mF:
            prev_idx = -1
            prev_line = None
            for line in rmF:
                idx = int(line.strip().split("_")[-2]) # account for the trailing '_temp'
                if idx == prev_idx:
                    print(f"Remove: {prev_line}")
                    meta_removed += 1            
                    prev_line = line
                    continue
                if prev_idx != -1:
                    if sarr[prev_idx-1] == 1:
                        mF.write(prev_line.strip()[:-5]+'\n') # get rid of the '_temp'
                        meta_kept += 1
                    else:
                        print(f"Remove: {prev_line}")
                        meta_removed += 1
                prev_line = line
                prev_idx = idx
            if sarr[prev_idx-1] == 1:
                mF.write(prev_line.strip()[:-5]+'\n') # get rid of the '_temp'
                meta_kept += 1
            else:
                print(f"Remove: {prev_line}")
                meta_removed += 1

    print(f"Meta lines kept/removed = {meta_kept}/{meta_removed}")

def clean_files(hdl, sarr, thr, sims_p_t):
    print(f"Cleaning tree files for {hdl} thread {thr}")

    begin = (thr-1)*sims_p_t+1 # inclusive
    end = thr*sims_p_t # inclusive
    file_kept = 0
    file_removed = 0
    no_file = 0
    for idx in range(begin, end+1):
        out_path = f"{hdl}/{hdl}_{idx}_samp.trees"
        if os.path.isfile(out_path):
            if sarr[idx-1] == 1:
                file_kept += 1
            else:
                file_removed += 1
                print(f"Remove: {out_path}")
                os.remove(out_path)
        else:
            print(f"File doesn't exist: {out_path}")
            no_file += 1
    print(f"Tree files kept/removed/none = {file_kept}/{file_removed}/{no_file}")

def main(args):
    if len(args) != 5:    #4 arguments
        return helpMsg

    mode = args[1]
    handle = args[2]
    no_sims_p_thr = int(args[3])
    no_thr = int(args[4])
    # indices are 1-based

    succ_arr = np.zeros(no_sims_p_thr*no_thr)
    #tree_arr = np.zeros(no_sims)

    with open(handle+".0exit", 'r') as eF:
        for line in eF:
            fld = line.strip().split('_')
            succ_arr[int(fld[0])-1] = 1
    print("Total success by EXITSTAT:", succ_arr.sum())

    if mode.isdigit():
        clean_files(handle, succ_arr, int(mode), no_sims_p_thr)
    else:
    	clean_meta_lines(mode, succ_arr)


    #print("Total tree files output:", tree_arr.sum())

    # mismatch = np.logical_xor(succ_arr, tree_arr)
    # print("Total mismatch:", mismatch.sum())

    # print("sim_no", "succ?", "tree?", sep='\t')
    # for idx in np.nonzero(mismatch)[0]:
    #     print(idx+1, succ_arr[idx], tree_arr[idx], sep='\t')        

    return 0


sys.exit(main(sys.argv))
