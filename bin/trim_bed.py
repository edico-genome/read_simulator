#!/usr/bin/env python2.7

import gzip, copy
import shutil
from collections import Counter

target_chrs='1'


bed = "/mnt/archive/varunj/MRJD_BedFiles/hs37d5.SegDupTable_sorted_v3.bed"
bed_out = "/mnt/archive/varunj/MRJD_BedFiles/hs37d5.SegDupTable_sorted_v3_chr1.bed"

print("Trim BED to small chromosome")
array_lines = []
array_regions = []


with open(bed) as stream:          
    for line in stream:
        # target_chrs can be an array or single string
        if line.split()[0] != target_chrs:
            continue

        # lines to keep
        array_lines.append(line)
        array_regions.append(line.split()[3])

# is this a vlrd bed? then we should only keep paired regions
hist = Counter(array_regions) # e.g. {"2": 3, "54": 3}
regions_to_keep = []
for key in hist:
    if hist[key] == 2:
        regions_to_keep.append(key)
regions_to_keep = set(regions_to_keep)

for i in reversed(range(0,len(array_lines))):
    if array_lines[i].split()[3] not in regions_to_keep:
        array_lines.pop(i)

if not array_lines:
    raise Exception("Empty target bed")


with open(bed_out, 'w') as stream:
    for line in array_lines:
        stream.write(line)
