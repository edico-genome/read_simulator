#!/usr/bin/env python2.7

import os

"""
@vlrd_normalnoise_errorrate_hg19_0.005_mason.000000001 contig=chr20 haplotype=0 length=101 orig_begin=57693380 orig_end=57693481 snps=0 indels=0 haplotype_infix=CAAGAGATTATTTTTTCTTCCTTAATATCCTCAGCGGATATTGTGCCTCAACGTGGCACAAAAATTTCTGAGTTGATGTCCAGATGAAATTTGCATGGAGA edit_string=MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM strand=forward
"""

base = "/mnt/archive/sim_data/dSim/vlrd_normalnoise_errorrate_hg19_0.005_mason/Mason"

file1_in = os.path.join(base, "vlrd_normalnoise_errorrate_hg19_0.005_mason_pi_0.001_pd_0.001_pmm_0.004_s_7451_N19032461_n_101_ll_410_le_22_1.fastq")
file2_in = os.path.join(base, "vlrd_normalnoise_errorrate_hg19_0.005_mason_pi_0.001_pd_0.001_pmm_0.004_s_7451_N19032461_n_101_ll_410_le_22_2.fastq")

file1_out=file1_in.replace(".fastq", "_cleaned.fastq")
file2_out=file2_in.replace(".fastq", "_cleaned.fastq")

print("\n\n")

for file_in, file_out in zip([file1_in, file2_in], [file1_out, file2_out]):
    print("{} -> {}".format(file_in, file_out))

    with open(file_in, 'r') as stream_in, open(file_out, 'w') as stream_out:
        idx = 0
        for line in stream_in:
            idx += 1
            line = line.replace("\n", "")

            # print(line)
            if line[0] == "@":
                line = line.split("snps")[0]
                line = line.replace("_normalnoise_errorrate_hg19_0.005_mason", "")
                line = line.replace("length=101", "")
                line = "_".join(line.split())

            if idx < 50:
                print line

            line = line + "\n"
            stream_out.write(line)

