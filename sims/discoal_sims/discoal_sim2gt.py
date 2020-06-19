#!/usr/bin/env python3

import os, sys
import subprocess
sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
 
import numpy as np
import utils
import pickle

helpMsg = '''
    usage: $./discoal_sim2gt.py <`swp` or `neu`> <handle> <s_low-s_high> <f_low-f_high> <c_low-c_high> <NOCHR> <thr> <sims_p_thr>

    Simulation de novo
    Run RELATE and output true and inferred tree sequence
    <s_low-s_high>: range of selection coefficient
    <f_low-f_high>: range of starting sweep frequency (0-0: hard sweep)
    <c_low-c_high>: range of ending sweep frequency

    Outputs:
        "discoal_pkl/discoal_<handle>_<NOCHR>_<`swp` or `neu`>_<thr>.pkl"
'''

## MACRO ##
DISCOAL_PATH = '/sonas-hs/siepel/hpc_norepl/home/mo/discoal-master'
MS_CMD = 'maiz_dem.txt'

## Simulation parameters ##
#Ne = 10000 # extract from ms command
pos = 0.5 #beneficial mutation positions (0-1)
length = 100000 #length of the simulated region
mu = 3e-8 # <= maize
#rho = 1.25e-8 # extract from ms command
SIGMA = 4 # time discretization for sweep simulation range: (4, 400),dt=(1/sigma*N)

def parse_msCmd(ms_cmd_ls, mu):

    theta = float(ms_cmd_ls[4])
    R = float(ms_cmd_ls[6])
    ms_win_size = int(ms_cmd_ls[7]) # i.e. nsites

    N0 = round(theta/(4*mu)/ms_win_size)
    rho = R/(4*N0)/ms_win_size

    dem_cmd = []
    for i in range(ms_cmd_ls.index("-eN"), len(ms_cmd_ls), 3):
        dem_cmd += ["-en", ms_cmd_ls[i+1], '0', ms_cmd_ls[i+2]]   

    return rho, N0, dem_cmd # list of population size change command for discoal

def parse_hdswp_cmd(cmd_line):
    
    args = cmd_line.split()
    N_curr = int(args[args.index('-N') + 1])
    alpha = float(args[args.index('-a') + 1])
    CAF = float(args[args.index('-c') + 1])
    SAF = float(args[args.index('-f') + 1])
    selcoef = alpha/(2*N_curr)
    
    return selcoef, SAF, CAF

def discoal2gt(discoal_file, cat, c_low, c_high):

    with open(discoal_file, "r") as discoalF:
        if cat == 'swp':
            SC, SAF, CAF = parse_hdswp_cmd(discoalF.readline().strip())
        elif cat == 'neu':
            SC = 0
            SAF = 0

        read_GT = False
        seek_onset = False # seek soft sweep onset
        onset_gen = -1
        for line in discoalF:
            if line[:8] == "segsites":
                segsites = int(line.strip().split()[1])
                gtm = np.empty((0, segsites), dtype=np.int8)
                continue
            if line[:9] == "positions":
                var_pos = np.array(line.strip().split()[1:], dtype=float)
                read_GT = True
                continue
            if line[:4] == "Freq":
                seek_onset = True
                continue
            if seek_onset:
                gbp_der_anc = line.strip().split()
                if len(gbp_der_anc) == 3 and float(gbp_der_anc[1]) < SAF: # 1st time point going backward
                    onset_gen = float(gbp_der_anc[0]) # in coalc. unit
                    seek_onset = False
                continue
            if read_GT:
                gtm = np.vstack((gtm, np.fromstring(line.strip(), dtype=np.int8) - ord("0")))

    gtm = np.transpose(gtm)
    if cat == 'neu':
        samp_idx = utils.samp_var(gtm, var_pos, c_low, c_high, 0.4, 0.6)
        foc_var_pos = var_pos[samp_idx]
        #foc_var_gt = gtm[samp_idx]
        CAF = np.mean(gtm[samp_idx])

    elif cat == 'swp':
        foc_var_pos = 0.5
        #foc_var_gt = gtm[var_pos == 0.5].flatten()

    #feaMtx = inf_fea(var_pos, gtm, foc_var_gt, foc_var_pos, str(2*10000), "2.5e-8", 1.25)
    
    return SC, SAF, CAF, onset_gen, foc_var_pos, var_pos, gtm

def main(args):
    if len(args) != 9:    #8 argument(s)
        return helpMsg
    
    mode = args[1] # <`swp` or `neu`>
    handle = args[2]
    s_low, s_high = list(map(float, args[3].split("-")))
    f_low, f_high = list(map(float, args[4].split("-")))
    c_low, c_high = list(map(float, args[5].split("-")))
    no_chrs = int(args[6])
    thr = int(args[7]) # one-based
    no_sims = int(args[8])

    with open(MS_CMD, "r") as inF:
        rho, Ne, discoal_dem = parse_msCmd(inF.readline().strip().split(), mu)

    theta = 4*Ne*mu*length
    R = 4*Ne*rho*length
    
    print(f"discoal_{handle}_{no_chrs}_{mode}_{thr}: {no_sims} sims", flush=True)

    SC_arr = np.exp(np.random.uniform(np.log(s_low), np.log(s_high), size=no_sims)) # selection coefficient
    SAF_arr = np.random.uniform(f_low, f_high, size=no_sims) # Starting AF
    CAF_arr = np.random.uniform(c_low, c_high, size=no_sims) # Current AF
    onset_arr = np.empty(no_sims)
    foc_var_pos_arr = np.empty(no_sims)
    var_pos_ls = []
    gtm_ls = []

    #wd = f"discoal_temp/discoal_{handle}_{no_chrs}_{mode}_{thr}"
    #temp_discoalF = f"discoal_{handle}_{no_chrs}_{mode}_{thr}_temp.discoal"

    for samp in range(no_sims):
        #if os.stat(discoalF_path).st_size == 0: continue
        temp_discoalF = f"discoal_temp/discoal_{handle}_{no_chrs}_{mode}_{thr}_{samp}.discoal"
        if mode == 'swp':
            sel = 2*Ne*SC_arr[samp]
            discoal_cmd = [f"{DISCOAL_PATH}/discoal", str(no_chrs), "1", str(length),
            "-t", str(theta), "-r", str(R),
            "-c", str(CAF_arr[samp]), "-f", str(SAF_arr[samp]), "-ws", "0", "-a", str(sel),
            "-N", str(Ne), "-i", str(SIGMA), "-x", str(pos), "-T"] + discoal_dem

        elif mode == 'neu':
            discoal_cmd = [f"{DISCOAL_PATH}/discoal", str(no_chrs), "1", str(length),
            "-t", str(theta), "-r", str(R), "-T"] + discoal_dem

        loop_cnt = 0
        while True:
            loop_cnt += 1
            print(f"SIM:{samp}, ATTEMPT:{loop_cnt}", flush=True)
            with open(temp_discoalF, "w") as outF:
                discoal_proc = subprocess.run(discoal_cmd, stdout=outF)
            if discoal_proc.returncode == 0: break

        SC, SAF, CAF, onset_gen, foc_var_pos, var_pos, gtm = discoal2gt(temp_discoalF, mode, c_low, c_high)
        
        if mode == 'neu':
            SC_arr[samp] = SC
            CAF_arr[samp] = CAF
            SAF_arr[samp] = SAF

        print(SC, f"{SAF}=>{CAF}", ";", SC_arr[samp], f"{SAF_arr[samp]}=>{CAF_arr[samp]}", ";", onset_gen) # sanity check
        onset_arr[samp] = onset_gen
        foc_var_pos_arr[samp] = foc_var_pos
        var_pos_ls.append(var_pos)
        gtm_ls.append(gtm)

        os.remove(temp_discoalF) # comment out for debugging

    # np.savez_compressed(f"discoal_meta/discoal_{handle}_{no_chrs}_{mode}_meta_{thr}",
    #     SC=SC_arr, SAF=SAF_arr, CAF=CAF_arr, onset=onset_arr, foc_var_pos=foc_var_pos_arr)
    with open(f"discoal_pkl/discoal_{handle}_{no_chrs}_{mode}_{thr}.pkl", "wb") as f:
        pickle.dump((SC_arr, SAF_arr, CAF_arr, onset_arr, foc_var_pos_arr, var_pos_ls, gtm_ls), f, pickle.HIGHEST_PROTOCOL)

    print("pickle fmt - (SC_arr, SAF_arr, CAF_arr, onset_arr, foc_var_pos_arr, var_pos_ls, gtm_ls)")
    print("Shape:", SC_arr.shape, SAF_arr.shape, CAF_arr.shape, onset_arr.shape, foc_var_pos_arr.shape, len(var_pos_ls), len(gtm_ls))

sys.exit(main(sys.argv))
