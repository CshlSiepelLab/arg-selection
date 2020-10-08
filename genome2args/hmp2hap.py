#!/usr/bin/env python3

import sys
import numpy as np
#import gzip

helpMsg = '''
        usage: $./hmp2hap.py <hmp_path> <out_hap_pref>

        Conversion from hapmap format (no ancestral info. for now)
'''

def parse_hmp_line(fields):
    
    #fields = hmp_line.strip().split()
    # 0:rs# 1:alleles 2:chrom 3:pos 4:strand 5:assembly# 6:center 7:protLSID 8:assayLSID 9:panelLSID 10:QCcode

    SNP_ID = fields[0]
    anc, der = fields[1].split("/") # assume ancestral allele is first
    chrm = int(fields[2])
    pos = int(fields[3])

    gt = [item for sublist in map(list, fields[11:]) for item in sublist]

    # mapping nucleotide to 0/1

    if gt.count(gt[0]) == len(gt): # monomorphic sites
        return None, None, 1, None

    match_al = np.sum(np.isin(gt, [anc, der]))
    if match_al == len(gt):
        code = 0
    else:
        code = 2
        print(f"{SNP_ID} m/m nucleotide: {len(gt)-match_al}", flush=True)

    nuc2bin = {anc:'0', der:'1'}
    gt_bin = list(map(lambda xx: nuc2bin.get(xx, '0'), gt)) # mismatch default to 0

    hap_line = " ".join([str(chrm), SNP_ID, str(pos), anc, der]+gt_bin)
    
    return SNP_ID, pos, code, hap_line

def main(args):
    if len(args) != 3:    #2 arguments
        return helpMsg

    hmp_path = args[1]
    hap_pref = args[2]

    with open(hmp_path, 'r') as hmpF:
        curr_chrm = -1
        cnt = [0, 0, 0] # bial, monomorphic, mismatch/missing

        for line in hmpF:
            if line[:3] == "rs#":
                indvls = line.strip().split()[11:]
                with open(hap_pref+".sample", "w") as sampF:
                    sampF.write("ID_1 ID_2 missing\n")
                    sampF.write("0 0 0\n")
                    for indvl in indvls:
                        sampF.write(f"{indvl} {indvl} 0\n")
            else:
                fields = line.strip().split()
                chrm = int(fields[2])
                if chrm != curr_chrm:
                    if curr_chrm > 0:
                        hapF.close()
                        metaF.close()
                        print(f"Chromosome {curr_chrm}: bial-{cnt[0]}, monomorph-{cnt[1]}, mis-{cnt[2]}", flush=True)
                        cnt = [0, 0, 0]

                    curr_chrm = chrm
                    hapF = open(f"{hap_pref}_{curr_chrm}.haps", "w")
                    metaF = open(f"{hap_pref}_{curr_chrm}.meta", "w")

                rs_ID, SNP_pos, code, hap_line = parse_hmp_line(fields)
                cnt[code] += 1

                if code == 1: continue # skip monomorphic sites
                hapF.write(hap_line+'\n')
                metaF.write(f"{rs_ID}\t{SNP_pos}\t{code}\n")
        hapF.close()
        metaF.close()
        print(f"Chromosome {curr_chrm}: bial-{cnt[0]}, monomorph-{cnt[1]}, mis-{cnt[2]}", flush=True)
    
    return 0

sys.exit(main(sys.argv))