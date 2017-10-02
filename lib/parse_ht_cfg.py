#!/usr/bin/env python

def parse_ht_config(cfg_file):

    d = {}
    nr_names = 0
    nr_len = 0

    with open(cfg_file) as stream:
        for line in stream:
            if line[0] == "#":
                continue
            if not line[0].split():
                continue

            firstCol = line.split()[0].strip()

            line = line.replace("\n", "")
            if "reference_sequences" == line.split()[0]:
                nr_sequences = int(line.split()[2])        
            if any((r == firstCol for r in
                    ["reference_len", "reference_len_raw", "reference_len_not_n"])):
                continue
            elif "reference_sequence" in line.split()[0]:
                name = line.split()[2]
                nr_names += 1
                print "name: ", name
            elif "reference_len" in line.split()[0]:
                print line
                r_len = line.split()[2]
                nr_len += 1
                d[name] = r_len
            if not nr_names == nr_len:
                raise Exception("parsing ht failed {} != {}".format(nr_names, nr_len))
    return d


