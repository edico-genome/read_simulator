#!/usr/bin/env python2.7
import subprocess
import os
import logging

logger = logging.getLogger(__name__)


##################################################
# CONSTANTS

script_dir = os.path.dirname(os.path.realpath(__file__))
rsv_script = os.path.join(script_dir, "unified_rsv_table.sh")
to_vcf_script = os.path.join(script_dir, "RSVSim_to_VCF.pl")
vcf_to_truth_table_script = os.path.join(script_dir, "rsvsim.vcf_to_cnv_truth_table.pl")

header = "##fileformat=VCFv4.2\n"
header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"

def run_cmd(cmd):
    subprocess.check_output(cmd, shell=True)

def get_order_file(fasta, outdir):
    fai_file = "{}.fai".format(fasta)
    assert os.path.isfile(fai_file)
    order_file = os.path.join(outdir, "order_file.txt")

    with open(fai_file) as stream_in, open(order_file, 'w') as stream_out:
        for line in stream_in:
            stream_out.write("{}\n".format(line.split()[0]))

    return order_file

def create_truth_files(fasta, outdir, csv_files):

    concatenated = os.path.join(outdir, "rsvsim.csv")
    with open(concatenated, 'w') as stream_out:
        for f in csv_files:
            f2 = os.path.join(outdir, f)
            with open(f2, 'r') as stream_in:
                for line in stream_in:
                    if line.replace("\n","").strip():
                        stream_out.write(line)

    # Order results
    order_file = get_order_file(fasta, outdir)
    
    normalized = os.path.join(outdir, "rsvsim_truth_output_normalized.txt")
    cmd = "{} < {} > {} ".format(rsv_script, concatenated, normalized)
    run_cmd(cmd)

    vcf = os.path.join(outdir, "RSVSim_truth.vcf")
    cmd = "echo {} > {}"
    cmd = cmd.format(header, vcf)
    run_cmd(cmd)

    cmd = "{} {} {} | sortBed -i /dev/stdin -faidx {} >> {}"
    cmd = cmd.format(to_vcf_script, normalized, fasta, order_file, vcf)
    run_cmd(cmd)
    logger.info("VCF Truth: {}".format(vcf))

    tsv = os.path.join(outdir, "RSVsim_truth.tsv")
    cmd = "{} {} > {}".format(vcf_to_truth_table_script, vcf, tsv)
    run_cmd(cmd)
    logger.info("TSV Truth: {}".format(tsv))

    return vcf, tsv
