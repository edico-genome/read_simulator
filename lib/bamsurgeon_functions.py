#!/usr/bin/env python2.7
import gzip

def vcf2bed(orig_vcf_path, new_vcf_path, bed_path, _logger):
    
    _logger.info("Create Bamsurgeon BED")
    _logger.info("Orig VCF: {}".format(orig_vcf_path))
    _logger.info("New VCF: {}".format(new_vcf_path))

    genotype = '0/1:1,1000:1000:10:10000,1000,0:0,1,1000,1000'
    
    with gzip.open(orig_vcf_path) as vcf_in, \
            open(bed_path, 'w') as bed_out, \
            open(new_vcf_path, 'w') as vcf_out: 
        idx = 0
        for line in vcf_in:

            # header
            if line[0] == '#':
                vcf_out.write(line)
                continue

            s = line.split()
            if any([line.split()[9] != "0/1", (len(s[4]) != 1)]):
                continue

            # new truth files
            vcf_out.write(line)
            bed_out.write("{}\t{}\t{}\t0.2\t{}\n".format(s[0],s[1],int(s[1])+1,s[4]))

            idx += 1

