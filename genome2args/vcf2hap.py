#!/usr/bin/env python3

import sys
import gzip

helpMsg = '''
        usage: $./vcf2hap.py <vcf_path> <out_hap_path>
'''

def parseVCFline(vcf_line):
    
    fields = vcf_line.strip().split()

    SNP_ID = fields[2]
    pos = int(fields[1])
    INFO = fields[7]
    gt = [item for sublist in map(lambda x: x.split("|"), fields[9:]) for item in sublist]
    if gt.count(gt[0]) == len(gt): # monomorphic sites
        return None, None, -3, None

    anc = fields[3]
    der = fields[4]

    if len(INFO) < 3:
        code = -1 # ancestral info n/a
    else:
        AA = INFO.split('=')[1][0]   
        if not AA.isalpha():
            code = -1 # ancestral info n/a
        
        elif AA.upper() == fields[3]: # REF is ancestral
            code = 0
            # anc = fields[3]
            # der = fields[4]
        elif AA.upper() == fields[4]: # ALT is ancestral
            code = 1
            anc = fields[4]
            der = fields[3]
            gt = list(map(lambda x: str(int(not int(x))), gt))
        else:
            code = -2 # mismatch to ancestral allele
    
    hap_line = " ".join([fields[0], fields[2], fields[1], anc, der]+gt)
    
    return SNP_ID, pos, code, hap_line

def main(args):
    if len(args) != 3:    #2 arguments
        return helpMsg

    vcf_path = args[1]
    hap_path = args[2]

    cnt = [0, 0, 0, 0, 0] # REF-anc, ALT-anc, monomorphic, mismatch, anc_n/a

    with gzip.open(vcf_path, 'rb') as vcfF:
        with open(hap_path, 'w') as hapF:
            with open(hap_path+".meta", 'w') as metaF:
                for line in vcfF:
                    string = line.decode()
                    if string[0] == '#': continue

                    rs_ID, SNP_pos, code, hap_line = parseVCFline(string)
                    cnt[code] += 1

                    if code == -3: continue # skip monomorphic sites
                    hapF.write(hap_line+'\n')
                    metaF.write(f"{rs_ID}\t{SNP_pos}\t{code}\n")

    print(f"Ref matching ancestral: {cnt[0]}\nAlt matching ancestral: {cnt[1]}\nMismatch: {cnt[-2]}\nAncestral n/a: {cnt[-1]}\nMonomorphic: {cnt[-3]}\n")
    
    return 0

sys.exit(main(sys.argv))