#!/usr/bin/env python2.7

"""
contig  start   stop    name
chr20   0       1000    target_1
chr20   1000    2000    target_2
chr20   2000    3000    target_3
chr20   3000    4000    target_4
"""

start = 0
stop = 249250621
chrom = "chr1"

start_region = 0
stop_region = 0

increment = 1000

print("contig  start   stop    name")
idx = 0
while ( stop_region < stop ):
    idx += 1
    name = "target_{}".format(idx)
    stop_region = min(stop_region + increment, stop)
    print "{}\t{}\t{}\t{}".format(chrom, start_region, stop_region, name)
    start_region += increment
