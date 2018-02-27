#!/usr/bin/env python2.7
import numpy as np
import collections

def get_dist():
    # dist = collections.OrderedDict()
    # "repeats"  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    # dist["1-1K"] = (1e3, [0.4, 0.15, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    # dist["1K-10K"] = (1e4, [0.45, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0,  0])
    # dist["10K-100K"] = (1e5, [0.6, 0.2, 0.1, 0.05, 0.05, 0, 0, 0, 0, 0])
    # dist["100K-Inf"] = (float('inf'), [0.7, 0.2, 0.1, 0, 0, 0, 0, 0, 0, 0])
    dist = [0.6, 0.4]
    return dist


def get_interval(x, line, dist):
    if x < 1:
        raise Exception("invalid event length in line: %s " % line)
    for key in dist:
        if x < dist[key][0]:
            return key

def resample_cnv_dup_repeats(file_in, file_out):
    dist = get_dist()
    arr = []
    with open(file_in) as stream_in, open(file_out, 'w') as stream_out:
        for idx, line in enumerate(stream_in):
            if idx == 0:
                stream_out.write(line)
                continue
            else:
                ss = line.split()
                try:
                    event_len = int(ss[4]) - int(ss[3])
                except:
                    print "unexpected line: ", line
                    raise

            #interval = get_interval(event_len, line, dist)
            # d = dist[interval][1]
            copies = str(np.random.choice(np.arange(1,3), p=dist))
            ss[6] = copies
            stream_out.write("\t".join(ss)+"\n")            
            print "len: %s %s" % (event_len, copies)

"""
if __name__ == "__main__":
    file_in = "/mnt/archive/sim_data/dSim/cnv_whg_sim_mixed_100bp_100kbp_20X_chr1"
    file_in += "/CNVgdbVCF/tandemDuplications.csv"
    file_out = "/staging/tmp/dups.csv"
    resample_cnv_dup_repeats(file_in, file_out)
    print file_in
    print file_out
"""

