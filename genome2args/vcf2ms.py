#!/usr/bin/env python3

import sys
import gzip
import numpy as np

helpMsg = '''
        usage: $./vcf2ms.py <vcf_path> <out_ms_path> <reg_start>
        Convert vcf file of a genomic loci to ms format
'''

N_HAPTYPE=182
REGION_LENGTH=2e6

def parseVCFline(vcf_line):
    
    fields = vcf_line.strip().split()
    INFO = fields[7]
    if len(INFO) < 3: return -1, None, None # ancestral info n/a
    AA = INFO.split('=')[1][0]
    gt = [int(item) for sublist in map(lambda x: x.split("|"), fields[9:]) for item in sublist]
    
    if not AA.isalpha(): return -1, None, None # ancestral info n/a
    
    if AA.upper() == fields[3]: # REF is ancestral
        code = 0
        # anc = fields[3]
        # der = fields[4]
    elif AA.upper() == fields[4]: # ALT is ancestral
        code = 1
        # anc = fields[4]
        # der = fields[3]
        gt = list(map(lambda x: int(not x), gt))
    else:
        return -2, None, None # mismatch to ancestral allele
    
    #hap_line = " ".join([fields[0], fields[2], fields[1], anc, der]+gt)
    
    return code, gt, int(fields[1]) # code, gt_ls, position

def main(args):
    if len(args) != 4:    #3 arguments
        return helpMsg

    vcf_path = args[1]
    ms_path = args[2]
    reg_start = int(args[3])

    geno_mtx = np.empty((0, N_HAPTYPE), dtype=np.int8)
    pos_ls = []
    cnt = [0, 0, 0, 0] # REF-anc, ALT-anc, mismatch, anc_n/a

    with gzip.open(vcf_path, 'rb') as vcfF:
        for line in vcfF:
            string = line.decode()
            if string[0] == '#': continue

            code, gt_ls, position = parseVCFline(string)
            cnt[code] += 1

            if code < 0: continue
            pos_ls.append((position - reg_start)/REGION_LENGTH)
            geno_mtx = np.vstack((geno_mtx, gt_ls))
            #hapF.write(hap_line+'\n')

    print(f"Ref matching ancestral: {cnt[0]}\nAlt matching ancestreal: {cnt[1]}\nMismatch: {cnt[-2]}\nAncestral n/a: {cnt[-1]}\n")
    print(f"Final gt mtx shape: {geno_mtx.shape}; pos ls length: {len(pos_ls)}", flush=True)

    geno_mtx = np.transpose(geno_mtx)

    with open(ms_path+".ms", 'w') as msF:
        msF.write(f"./discoal {N_HAPTYPE} 1 {int(REGION_LENGTH)}\n1 1\n\n//\n")
        msF.write(f"segsites: {len(pos_ls)}\n")
        print("positions:", *pos_ls, file=msF)
        for j in range(0, np.shape(geno_mtx)[0]):
            print(*geno_mtx[j], sep="", file=msF)

    return 0

sys.exit(main(sys.argv))