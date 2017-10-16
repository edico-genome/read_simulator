#!/usr/bin/env python 

target_bed = "/mnt/archive/gavinp/1000_genomes/20120518.consensus_add50bp.chrom.bed"
target_bed_contigs_name = target_bed.replace(".bed", "_contigs_named.bed")

if __name__ == "__main__":
    with open(target_bed) as stream_in, open(target_bed_contigs_name, 'w') as stream_out:
        stream_out.write("contig\tstart\tstop\tname\n")
        for idx, line in enumerate(stream_in):
            new_line = "{}\ttarget-{}\n".format(line.replace("\n",""), idx+1)
            stream_out.write(new_line)

print("Wrote results to: {}".format(target_bed_contigs_name))
