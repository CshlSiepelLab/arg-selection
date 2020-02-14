#!/usr/bin/env python3

import sys

helpMsg = '''
        usage: $./make_poplabels.py <samp_file> <pop_label> <group_label> <out_pref>
'''

def main(args):
	if len(args) != 5:    #4 arguments
		return helpMsg

	samp_file=args[1]
	pop_label=args[2]
	group_label=args[3]
	out_pref=args[4]
	
	with open(samp_file, 'r') as sampF:
		with open(out_pref+".poplabels", 'w') as popF:
			popF.write("sample population group sex\n")
			
			sampF.readline()
			sampF.readline()
			for line in sampF:
				fields = line.strip().split()
				popF.write(" ".join((fields[0], pop_label, group_label, '1'))+"\n")

	return 0

sys.exit(main(sys.argv))