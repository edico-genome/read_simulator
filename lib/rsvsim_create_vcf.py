#!/usr/bin/env python2.7
import subprocess
import os

##################################################
# CONSTANTS

script_dir = os.path.dirname(os.path.realpath(__file__))
rsv_script = os.path.join(script_dir, "../bin", "unified_rsv_table.sh")
to_vcf_script = os.path.join(script_dir, "../bin", "RSVSim_to_VCF.pl")


def run_cmd(cmd, _logger):
    _logger.info(cmd)
    subprocess.check_output(cmd, shell=True)

def get_order_file(fasta, outdir, _logger):
    fai_file = "{}.fai".format(fasta)
    if not os.path.isfile(fai_file):
        cmd = "samtools faidx {}".format(fasta)
        run_cmd(cmd, _logger)

    order_file = os.path.join(outdir, "order_file.txt")
    with open(fai_file) as stream_in, open(order_file, 'w') as stream_out:
        for line in stream_in:
            stream_out.write("{}\n".format(line.split()[0]))

    return order_file

def sort_tsv(tsv, order_file, tsv_sorted, _logger):

    # write header
    with open(tsv) as stream_in, open(tsv_sorted, 'w') as stream_out:
        for line in stream_in:
            stream_out.write(line)
            break

    # append sorted tsv body
    # have to print col 2 twice to get sortBed to think this is a bed file
    cmd = "awk '{ print $1\"\\t\"$2\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6 }' " + tsv
    cmd += "| sortBed -i /dev/stdin -faidx " + order_file + " | cut -f 1-2,4-7 >> " + tsv_sorted
    _logger.info(cmd)
    
    try:
        res = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        _logger.info("TSV Truth: {}".format(tsv_sorted))
    except:
        _logger.error("Failed to sort tsv")
        raise

def create_truth_files(fasta, outdir, csv_files, _logger):

    concatenated = os.path.join(outdir, "rsvsim.csv")
    with open(concatenated, 'w') as stream_out:
        for f in csv_files:
            f2 = os.path.join(outdir, f)
            with open(f2, 'r') as stream_in:
                for line in stream_in:
                    if line.replace("\n","").strip():
                        stream_out.write(line)

    # Order results
    order_file = get_order_file(fasta, outdir, _logger)
    
    normalized = os.path.join(outdir, "rsvsim_truth_output_normalized.txt")
    cmd = "{} < {} > {} ".format(rsv_script, concatenated, normalized)
    run_cmd(cmd, _logger)

    """
    #CHROM  START   END     ID      CHROM_2 START_2 END_2   LENGTH  REPEATS
    chr20   147379  182086  deletion1       NA      NA      NA      34708   1
    chr20   622800  623205  deletion2       NA      NA      NA      406     1
    chr20   939636  939822  deletion3  
    """

    # remove any remaining overlaps
    var_dict = {}
    idx = 0
    with open(normalized, 'r') as stream_in:
        for line in stream_in:
            idx += 1
            if idx == 1:
                header_line = line
            else:
                chrom = line.split()[0]
                start = int(line.split()[1])
                if chrom not in var_dict:
                    var_dict[chrom] = {}
                var_dict[chrom][start] = line

    assert (len(var_dict.keys()) == 1), "only support 1 target chr at present"
    first_chr = var_dict.keys()[0]
    start_pos = sorted(var_dict[first_chr])

    normalized_overlaps_removed = os.path.join(outdir, "rsvsim_truth_output_normalized_overlaps_removed.txt")
    with open(normalized_overlaps_removed, 'w') as stream_out:
        stream_out.write(header_line)
        idx = 0
        for sp in start_pos:
            idx += 1
            line = var_dict[first_chr][sp]
            if idx == 1:
                stream_out.write(line)
                _prev_var_end = int(line.split()[2])
            else:
                _var_start = int(line.split()[1])
                if _var_start > _prev_var_end:
                    stream_out.write(line)
                    _prev_var_end = int(line.split()[2])
                    _logger.info('prev var end {}, this var start: {}'.format(_prev_var_end, _var_start))
                else:
                    _logger.info('discarding overlapping variant: {}'.format(line))

                
    vcf = os.path.join(outdir, "RSVSim_truth.vcf")
    vcf_gz = os.path.join(outdir, "RSVSim_truth.vcf.gz")
    with open(vcf, 'w') as stream:
        stream.write("##fileformat=VCFv4.2\n")
        stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

    cmd = "{} {} {} | sortBed -i /dev/stdin -faidx {} >> {}"
    cmd = cmd.format(to_vcf_script, normalized_overlaps_removed, fasta, order_file, vcf)
    run_cmd(cmd, _logger)
    _logger.info("VCF Truth: {}".format(vcf))

    # create truth tsv
    vcf_to_truth_table_script = os.path.join(
        script_dir, "../bin",
        "rsvsim.vcf_to_cnv_truth_table.pl")

    tsv = os.path.join(outdir, "RSVsim_truth.tsv")
    cmd = "{} {} > {}".format(vcf_to_truth_table_script, vcf, tsv)
    run_cmd(cmd, _logger)
    
    # sort tsv
    tsv_sorted = os.path.join(outdir, "RSVsim_truth_sorted.tsv")
    sort_tsv(tsv, order_file, tsv_sorted, _logger)
        
    # compress vcf for bcftools
    cmd = "bgzip -c {} > {}".format(vcf, vcf_gz)
    run_cmd(cmd, _logger)

    cmd = "tabix {}".format(vcf_gz)
    run_cmd(cmd, _logger)

    return vcf_gz, tsv_sorted
