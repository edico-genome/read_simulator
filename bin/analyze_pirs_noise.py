#!/usr/bin/env python2.7

import sys
import os 
import gzip

for length, f in ((100, "/mnt/archive/sim_data/dSim/dsim_masonNoise_10X_chr1_100bp_4xIns/Pirs/pirs_100_400.read.info.gz"),
                  (400, "/mnt/archive/sim_data/dSim/dsim_masonNoise_10X_chr1_400bp_4xIns/Pirs/pirs_400_1600.read.info.gz")):


    print '=================== {} ========================='.format(f.split('/')[-1])

    print f
    print length

    counts = [0]*length
    nr_lines = 0

    with gzip.open(f) as stream:
        for line in stream:

            if line[0] == "#":
                continue

            line = line.replace("\n", "")
            nr_lines += 1

            if nr_lines > 5000000:
                count_ratio = [ i*1.0/nr_lines for i in counts ]
                count_ratio_cut = [ "%.4f" % i for i in count_ratio ]
                print count_ratio_cut
                break

            mismatches = line.split()[7]
            if mismatches == "-":
                continue
            mismatches = mismatches.split(":")

            for m in mismatches:
                idx = int(m.split(",")[0])
                counts[idx-1] += 1

