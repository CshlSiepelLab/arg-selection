#!/usr/bin/env python3

import sys
import gzip

helpMsg = '''
        usage: $./vcf2hap.py <vcf_path> <out_hap_path>
'''

def parseVCFline(vcf_line):
    
    fields = vcf_line.strip().split()
    INFO = fields[7]
    if len(INFO) < 3: return -1, None # ancestral info n/a
    AA = INFO.split('=')[1][0]
    gt = [item for sublist in map(lambda x: x.split("|"), fields[9:]) for item in sublist]
    
    if not AA.isalpha(): return -1, None # ancestral info n/a
    
    if AA.upper() == fields[3]: # REF is ancestral
        code = 0
        anc = fields[3]
        der = fields[4]
    elif AA.upper() == fields[4]: # ALT is ancestral
        code = 1
        anc = fields[4]
        der = fields[3]
        gt = list(map(lambda x: str(int(not int(x))), gt))
    else:
        return -2, None # mismatch to ancestral allele
    
    hap_line = " ".join([fields[0], fields[2], fields[1], anc, der]+gt)
    
    return code, hap_line

def main(args):
    if len(args) != 3:    #2 arguments
        return helpMsg

    vcf_path = args[1]
    hap_path = args[2]

    cnt = [0, 0, 0, 0] # REF-anc, ALT-anc, mismatch, anc_n/a

    with gzip.open(vcf_path, 'rb') as vcfF:
        with open(hap_path, 'w') as hapF:
            for line in vcfF:
                string = line.decode()
                if string[0] == '#': continue

                code, hap_line = parseVCFline(string)
                cnt[code] += 1

                if code < 0: continue
                hapF.write(hap_line+'\n')

    print(f"Ref matching ancestral: {cnt[0]}\nAlt matching ancestreal: {cnt[1]}\nMismatch: {cnt[-2]}\nAncestral n/a: {cnt[-1]}\n")
    
    return 0

sys.exit(main(sys.argv))