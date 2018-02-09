#!/usr/bin/env python2.7

f = "/mnt/archive/varunj/MRJD_BedFiles/hs37d5.SegDupTable_sorted_v3_chr1.bed"

counter = 0
diff_total = 0


with open(f) as stream:
    for line in stream:
        if not line:
            continue
        counter += 1
        _from = int(line.split()[1])
        _to = int(line.split()[2])
        diff = _to - _from 
        diff_total += diff

print "lines: ", counter
print "avg dist: ", (diff_total / counter)
