#! /usr/bin/env python

import sys, os
import random

helpMsg = '''
usage: $./samp_sc_af.py sc_low sc_high af_low af_high replPerThr noThr sim_TAG
'''

def main(args):
	if len(args) != 8:    #7 arguments
		return helpMsg

	minSC = float(args[1])
	maxSC = float(args[2])
	minAF = float(args[3])
	maxAF = float(args[4])

	replPerThr = int(args[5])
	noThr = int(args[6])
	TAG = args[7]

	#os.mkdir("discoal_"+TAG)

	simID_0b = 0

	for thr in range(1, noThr+1):
		with open("discoal_"+TAG+"_"+str(thr)+"_id_sc_af.txt", "w") as outF:
			for sim in range(replPerThr):

				ID_sc_af = [str(simID_0b), str(round(random.uniform(minSC,maxSC), 4)), str(round(random.uniform(minAF,maxAF), 2))]
				outF.write("\t".join(ID_sc_af)+"\n")
				simID_0b += 1

	return 0

sys.exit(main(sys.argv))