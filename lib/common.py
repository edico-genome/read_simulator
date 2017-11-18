import logging
import pdb, sys, os
import subprocess

logger = logging.getLogger(__name__)

###########################################################
class MPdb(pdb.Pdb):
    """A Pdb subclass that may be used
    from a forked multiprocessing child

    """
    def interaction(self, *args, **kwargs):
        _stdin = sys.stdin
        try:
            sys.stdin = open('/dev/stdin')
            pdb.Pdb.interaction(self, *args, **kwargs)
        finally:
            sys.stdin = _stdin


###########################################################
def run_process(cmd, _logger, outfile=None):
    """ helper function to run cmd """

    if not isinstance(cmd, list):
        cmd = cmd.split()
    arguments = cmd

    _logger.info("Run cmd: {}".format(" ".join(arguments)))

    if outfile:
        with open(outfile, 'w') as f:
            process = subprocess.Popen(
                arguments,
                stdout=f)
            process.wait()
            _logger.info("wrote output to: {}".format(outfile))
    else:
        process = subprocess.Popen(
            arguments,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    process.wait()
    output, error = process.communicate()
    _logger.info(output)
    
    _logger.info("cmd exit code: {}".format(process.returncode))
    if process.returncode:
        _logger.error("cmd: {}".format(" ".join(arguments)))
        raise PipelineExc("bash cmd failed: {}".format(error))

    return output


###########################################################
# exceptions

class PipelineExc(Exception):
    """
    only a specific pipeline failed, still continue with next pipelines
    """
    pass


class Fatal(Exception):
    """
    exception in the dsim run, all pipelines abort
    """
    pass


###########################################################
def trim_fasta(_module):
    """
    if we should only use a subset of the input fasta
    useful e.g. if we want to speed up analyses
    """
    target_chrs = _module.module_settings['target_chrs']

    _module.logger.info("Trim Fasta for chromosome(s): {}".format(target_chrs))
    # workdir - no need to keep this fasta
    out_fasta = os.path.join(
        _module.module_settings['workdir'],
        "chr_{}.fa".format(target_chrs))

    # trim the fasta
    cmd = ["samtools", "faidx",
           _module.module_settings['fasta_file'],
           target_chrs]
    
    run_process(cmd, _module.logger, out_fasta)

    # index the new fasta
    cmd = ["samtools", "faidx", out_fasta]
    run_process(cmd, _module.logger)

    # update db and pipeline settings
    _module.module_settings["fasta_file"] = out_fasta


def trim_vcf(_module):
    """
    if we should only use a subset of the input vcf
    useful e.g. if we want to speed up analyses
    """
    target_chrs = _module.module_settings['target_chrs']

    _module.logger.info("Trim VCF for chromosome(s): {}".format(target_chrs))
    out_vcf = os.path.join(
        _module.module_settings['outdir'],
        "chr_{}.vcf.gz".format(target_chrs))

    cmd = ["bcftools", "filter", "--output-type", "z",
           "--regions", target_chrs,
           _module.module_settings['truth_set_vcf']]
    run_process(cmd, _module.logger, out_vcf)

    cmd = ["bcftools", "index", out_vcf]
    run_process(cmd, _module.logger)

    _module.module_settings["truth_set_vcf"] = out_vcf


###########################################################
def trim_bam(_module):
    """
    if we should only use a subset of the input vcf
    useful e.g. if we want to speed up analyses
    """
    in_bam = _module.module_settings.get("dragen_BAM", None)
    if not in_bam:
        return 

    target_chrs = _module.module_settings['target_chrs']
    out_bam = os.path.join(
        _module.module_settings['outdir'],
        "chr_{}.bam".format(target_chrs))

    _module.logger.info("Trim BAM for chromosome(s): {}".format(target_chrs))

    cmd = ["samtools", "view", "-b", "-h", in_bam, target_chrs]
    run_process(cmd, _module.logger, out_bam)

    cmd = ["samtools", "index", out_bam]
    run_process(cmd, _module.logger)

    _module.module_settings["dragen_BAM"] = out_bam
    _module.db_api.upload_to_db('dragen_BAM', out_bam)
  

###########################################################
def parse_ht_config(cfg_file):
    """
    a function to parse the ht config file and extract sequence length
    """
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


###########################################################
def remove_contig_name_descriptions(_module):
    """ some tools (e.g. RSVSim/ bcftools) fail when we use fastas
    with long contig names like 
    >1 dna:chromosome chromosome:GRCh37:1:1:249250621:1 
    this function strips the crud and copies the file to staging """

    # only have one process create this file
    _module.lock.acquire()
    basename = os.path.basename(_module.module_settings["fasta_file"])
    new_fasta = os.path.join(_module.module_settings["shared_dir"], "{}_mod".format(basename))
    
    # look at first line to determine if we need to process
    with open(_module.module_settings["fasta_file"], 'r') as stream_in:
        first_line = stream_in.readline().strip("\n")
        if len(first_line.split()) == 1:
            _module.logger.info("Fasta file in expected format, no need to process")
            _module.lock.release()
            return
        
    _module.logger.info("Fasta file first line: {}".format(first_line))
    _module.logger.info("Fasta file not in expected format")

    if os.path.isfile(new_fasta):
        _module.logger.info("Fasta previously preprocessed and available for use by RSVSim")
    else:
        _module.logger.info("Preprocessing fasta file for use by RSVSim")            
        with open(_module.module_settings["fasta_file"], 'r') as stream_in, \
             open(new_fasta, "w") as stream_out:
            try:
                for line in stream_in:
                    if line[0] == ">":
                        line = "{}\n".format(line.split()[0])
                    stream_out.write(line)
            except:
                _module.logger.error(
                    "Failed to preprocess fasta line: {}".format(line))
                raise
                    
        cmd = "samtools faidx {}".format(new_fasta)
        _module.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    _module.lock.release()
    _module.module_settings["fasta_file"] = new_fasta


###########################################################
def sort_target_region_bed(_module):
    sorted_bed = os.path.join(_module.module_settings["outdir"],
                              "target_sorted.bed")
    cmd = "sort -V -k1,1 -k2,2g {} | grep -v GL > {}".format(
        _module.module_settings["target_region_bed"], sorted_bed)
    
    try:
        _module.logger.info(cmd)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    except Exception as e:
        raise PipelineExc(e)
        
    _module.module_settings["sorted_bed"] = sorted_bed
