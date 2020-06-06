#! /usr/bin/env python

import sys, os
import random

helpMsg = '''
usage: $./samp_discrete_sc_af.py replPerThr noThr sim_TAG
'''

def main(args):
	if len(args) != 4:    #3 arguments
		return helpMsg

	## Hard-coded ##
	SC_range = [0.001, 0.0025, 0.005, 0.0075, 0.01]
	minAF = 0.1
	maxAF = 0.95

	replPerThr = int(args[1])
	noThr = int(args[2])
	TAG = args[3]

	#os.mkdir("discoal_"+TAG)

	simID_0b = 0

	for thr in range(1, noThr+1):
		with open("discoal_"+TAG+"_"+str(thr)+"_id_sc_af.txt", "w") as outF:
			for sim in range(replPerThr):

				ID_sc_af = [str(simID_0b), str(random.choice(SC_range)), str(round(random.uniform(minAF,maxAF), 2))]
				outF.write("\t".join(ID_sc_af)+"\n")
				simID_0b += 1

	return 0

sys.exit(main(sys.argv))