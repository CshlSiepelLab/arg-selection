#!/usr/bin/env python3

import os, sys
import subprocess
import numpy as np
#import pandas as pd
import allel
import tskit

RELATE_PATH = '/sonas-hs/siepel/hpc_norepl/home/mo/relate_v1.0.17_x86_64_static/'
#RELATE_PATH = '~/relate_v1.0.16_MacOSX/'
ii32MAX = np.iinfo(np.int32).max

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
            print("Milestone: attempt exceeds 20, ABORTS!!", file=sys.stderr)
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

def calc_H1(gt_mtx):

    #pos = np.around(pos * 100000) # convert position to coordinate in 100kb region
    hArr = allel.HaplotypeArray(gt_mtx)
    acArr = hArr.count_alleles() 

    H1, H12, H123, H2H1 = allel.garud_h(hArr)

    return H1
