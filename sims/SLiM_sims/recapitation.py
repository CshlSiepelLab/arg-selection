#!/usr/bin/env python3

import sys, os
import msprime
import pyslim
import numpy as np

helpMsg = '''
        usage: $./recapitation.py <slim_tree.trees> <slim_dem.params> <N_chrs> <out_pref>

        - recapitates and samples from the tree-sequence output from SLiM simulations
'''

def main(args):
    if len(args) != 5:    #4 arguments
        return helpMsg

    slim_tree_path = args[1]
    slim_params_path = args[2]
    N = int(args[3])
    out_path = args[4]+".trees"

    with open(slim_params_path, "r") as paramF:
        lines = paramF.readlines()

    mu = float(lines[0].strip())
    rho = float(lines[1].strip())
    Ne_recap = int(lines[2].strip().split()[1])

    ts = pyslim.load(slim_tree_path)
    ts_recap = ts.recapitate(recombination_rate=rho, Ne=Ne_recap)
    sampN = np.random.choice(ts_recap.samples(), size=N, replace=False)
    ts_samp = ts_recap.simplify(samples=sampN)

    ts_samp.dump(out_path)

    return 0

sys.exit(main(sys.argv))