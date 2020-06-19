#!/usr/bin/env python3

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../..")) # directory of module fea_util
from contextlib import contextmanager

import msprime, pyslim
import numpy as np
import tszip

from RELATE_util import run_RELATE

helpMsg = '''
        usage: $./slimSims2ts.py <.param file> <scale> <max_no_sims> <thr> <tot_thr> <inPref> <outPref>
            Take slim simulation output, run RELATE and save both true and inferred tree sequence file for each sim
            <inPref> and <outPref> should be directory names!
'''

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

def sim2ts(ts_mut, mu, rho_cMpMb, N0):

    # ts_samp = pyslim.load(sim_path)
    # ts_mut = msprime.mutate(ts_samp, rate=scaled_mu, keep=True) # random_seed=958,

    ppos_ls = get_site_ppos(ts_mut)
    GTM = ts_mut.genotype_matrix()
    
    mask = (ppos_ls != -1)
    ppos_ls = ppos_ls[mask]
    GTM = GTM[mask, :]

    ts_inf, _ = run_RELATE(ppos_ls, GTM, str(2*N0), rho_cMpMb=rho_cMpMb, mut_rate=str(mu))

    return ts_inf # inferred tree-seqs

# @contextmanager is just an easier way of saying cd = contextmanager(cd)
@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def main(args):
    if len(args) != 8:    #7 arguments
        return helpMsg

    parent_dir = os.getcwd() # directory where the job is submitted
    param_path = args[1]
    scale = int(args[2])
    no_sims = int(args[3])
    thr = int(args[4]) # should be 1-indexed
    tot_thr = int(args[5])
    inPref = args[6]
    outPref = args[7]

    with open(param_path, "r") as paramF:
        lines = paramF.readlines()
    scaled_mu = float(lines[0].strip())
    scaled_rho = float(lines[1].strip())
    scaled_N0 = int(lines[-1].strip().split()[1])

    mu = scaled_mu/scale
    rho_cMpMb = scaled_rho/scale*100*1e6 # 1cM = 1e-2 crossover
    N0 = scaled_N0*scale

    # if mode == 's':
    #     meta_data = np.genfromtxt(metaF, usecols=(1, 2, 3, 4), dtype=None)
    #     # e.g. `%% 0.100992 2347 0.520421 SLiM_trial_swp/SLiM_trial_swp_4501`
    #     no_sims = meta_data.shape[0]

    tasks = no_sims//tot_thr
    a_idx = (thr-1)*tasks # inclusive
    if thr == tot_thr:
        b_idx = no_sims # exclusive
    else:
        b_idx = thr*tasks # exclusive
    # indices are 0-based

    print("Processing: [", a_idx, b_idx, ")", flush=True)

    wd = outPref+'/'+'RELATE_temp_'+str(thr)
    os.mkdir(wd, 0o755)

    # idx_ls = []
    # sc_ls = []
    # onset_ls = []
    # caf_ls = []
    #fv_idx_ls = []
    #cnt = 0

    log_f = open(outPref+"_"+str(thr)+".log", 'a')
    with cd(wd):
        for r_idx in range(a_idx, b_idx):
            # if mode == 'n':
            ID = r_idx + 1 # convert 0-based index to 1-based index
            sim_path = parent_dir+"/"+inPref+"/"+inPref+"_"+str(ID)+"_samp.trees"
            # elif mode == 's':
            #     ID = int(meta_data[r_idx][3].split(b'_')[-1]) # retrieve 1-based index from meta file
            #     sim_path = parent_dir+"/"+meta_data[r_idx][3].decode()+"_samp.trees"

            if not os.path.isfile(sim_path):
                print(ID, "SKIPPED", file=log_f, flush=True)
                continue
            # if mode == 's':
            #     idx_ls.append(ID)
            #     sc_ls.append(meta_data[r_idx][0])
            #     onset_ls.append(meta_data[r_idx][1])
            #     caf_ls.append(meta_data[r_idx][2])

            outFP = parent_dir+"/"+outPref+"/"+outPref+"_"+str(ID)
            if os.path.isfile(outFP+"_tru.trees") and os.path.isfile(outFP+"_inf.trees.tsz"): continue

            print("Input:", sim_path, flush=True)
            ts_samp = pyslim.load(sim_path)
            ts_tru = msprime.mutate(ts_samp, rate=scaled_mu, keep=True)
            ts_inf = sim2ts(ts_tru, mu, rho_cMpMb, N0)

            ts_tru.dump(outFP+"_tru.trees")
            #ts_inf.dump(outFP+"_inf.trees")
            tszip.compress(ts_inf, outFP+"_inf.trees.tsz")
            print(ID, "SUCCESS", file=log_f, flush=True)

    # if mode == 's': np.savez_compressed(parent_dir+"/"+outPref+"/"+outPref+"_meta_"+str(thr), 
    #     idx=idx_ls, sc=sc_ls, onset=onset_ls, caf=caf_ls)
    log_f.close()
    os.rmdir(wd)

    return 0

sys.exit(main(sys.argv))
