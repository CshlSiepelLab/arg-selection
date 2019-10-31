#!/usr/bin/env python3

import sys
import gzip

helpMsg = '''
        usage: $./make_uniform_map.py <hap_file> <rrate_cM/Mb> <out_pref>
'''

def main(args):
	if len(args) != 4:    #3 arguments
		return helpMsg

	hap_file=args[1]
	rrate=float(args[2])
	out_pref=args[3]
	
	#with gzip.open(hap_file, 'rb') as hapF:
	with open(hap_file, 'r') as hapF:
		with open(out_pref+".map", 'w') as mapF:
			ppos = 0
			rdist = 0
			mapF.write("pos COMBINED_rate Genetic_Map\n")
			mapF.write(" ".join((str(ppos), str(rrate), str(rdist)))+"\n")
			
			for line in hapF:
				fields = line.strip().split()
				next_ppos = int(fields[2])
				rdist = rdist + rrate/1e6*(next_ppos - ppos)
				ppos = next_ppos
				mapF.write(" ".join((str(ppos), str(rrate), str(rdist)))+"\n")

	return 0

sys.exit(main(sys.argv))
