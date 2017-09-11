"""
A module defines a unit of work ( creating a truth VCF/ fasta etc )
"""

import os
import sys
import glob
import logging
from vlrd import vlrd_functions
import subprocess
from modules_base import ModuleBase

logger = logging.getLogger(__name__)
this_dir_path = os.path.dirname(os.path.realpath(__file__))


###########################################################
class VLRDVCF(ModuleBase):
    default_settings = {}
    expected_settings = ["varrate"]

    def run(self):
        # get last minute settings
        self.get_dataset_ref()
        self.get_target_bed()

        for key in ["ref_type", "ref_fasta", "target_bed"]:
            if not self.module_settings[key]:
                print "VLRD, dataset {} missing: {}".format(self.dataset_name, key)
                sys.exit(1)

        # run
        res = vlrd_functions.create_truth_vcf_and_fastas(self.module_settings)
        logger.info(res)
        self.db_api.set_fastas(res['fasta0'], res['fasta1'])
        self.db_api.post_truth_vcf(res['truth_vcf'])


###########################################################
class AltContigVCF(ModuleBase):
    default_settings = {}
    expected_settings = ["contig", "alt-contig-1", "alt-contig-2", "alt-sam"]

    def get_alt_contig_lengths(self):
        """ lookup alt contig length in dict - fix this mock """
        logger.info("Get alt contig lengths")

        with open(self.module_settings["dict_file"], 'r') as stream:
            for line in stream:
                len_col = line.split()[2]
                if "LN" in len_col:
                    for c_idx in ["alt-contig-1", "alt-contig-2"]:
                        contig = self.module_settings[c_idx]
                        if contig in line:
                            try:
                                c_len = int(len_col.replace("LN:", ""))
                                self.module_settings[c_idx + "-len"] = c_len
                                logger.info("%s: %s" % (self.module_settings[c_idx], self.module_settings[c_idx+"-len"]))
                            except Exception as e:
                                logger.error("Contig len not found/ integer in line: {}".format(line))

        for i in ["alt-contig-1-len", "alt-contig-2-len"]:
            assert isinstance(self.module_settings[i], int), "Invalid contig len detected from dict"

    def convert_alt_to_pri_coords(self):
        """ get coordinates which we'll use for the truth vcf generation """
        logger.info("Convert the alternate contig to primary coords")
        alt_to_pri_script = os.path.join(this_dir_path, "alt_contig", "convert_alt_to_pri_coords.pl")
        alt_contig_sam = "/opt/bwakit-0.7.12-0/resource-GRCh38/hs38DH.fa.alt"

        cmd = "{} {} {} 1 {} 2> /dev/null".format(alt_to_pri_script,
            alt_contig_sam,
            self.module_settings["alt-contig-1"],
            self.module_settings["alt-contig-1-len"])
        logger.info("Run cmd: {}".format(cmd))
        res = subprocess.check_output(cmd, shell=True)
        logger.info("Response: {}".format(res))
        try:
            _contig, _from, _to = res.split()[0:3]
            self.module_settings["main_contig_from"] = _from
            self.module_settings["main_contig_to"] = _to
        except Exception as e:
            logger.error('Failed to align alt-contig-1 to main contig')
            sys.exit(1)

    def create_modified_fastas(self):
        fa0 = "{outdir}/mod_fasta_0.fa".format(**self.module_settings)
        fa1 = "{outdir}/mod_fasta_1.fa".format(**self.module_settings)

        cmd = "samtools faidx {fasta_file} {alt-contig-1}:1-{alt-contig-1-len} > %s" % fa0
        cmd = cmd.format(**self.module_settings)
        subprocess.check_call(cmd, shell=True)
        logger.info("create modified fasta 0: cmd: {}".format(cmd))

        cmd = "samtools faidx {fasta_file} {alt-contig-2}:1-{alt-contig-2-len} > %s" % fa1
        cmd = cmd.format(**self.module_settings)
        subprocess.check_call(cmd, shell=True)
        logger.info("create modified fasta 1: cmd: {}".format(cmd))
        self.db_api.set_fastas(fa0, fa1)

    def fasta_sam_to_vcf(self):
        script_path = os.path.join(this_dir_path, "alt_contig", "fasta_sam_to_vcf.pl")
        self.module_settings["truth_vcf"] = os.path.join(self.module_settings["outdir"], "alt_contig_truth.vcf")
        cmd = script_path + " {fasta_file} {alt-sam} {alt-contig-1}:1-{alt-contig-1-len} {alt-contig-2}:1-{alt-contig-2-len} > {truth_vcf}"
        cmd = cmd.format(**self.module_settings)
        logger.info("create truth vcf cmd: {}".format(cmd))
        try:
            subprocess.check_call(cmd, shell=True)
        except Exception as e:
            raise Exception('Failed to create truth VCF: {}'.format(e))
        self.db_api.post_truth_vcf(self.module_settings["truth_vcf"])

    def run(self):
        self.before_run()
        self.get_alt_contig_lengths()
        self.convert_alt_to_pri_coords()
        self.create_modified_fastas()
        self.fasta_sam_to_vcf()
        logger.info(" - done")


###########################################################
class Pirs(ModuleBase):
    default_settings = {
        "PE100": "/opt/pirs-2.0.1/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz",
        "indels": "/opt/pirs-2.0.1/Profiles/InDel_Profiles/phixv2.InDel.matrix",
        "gcdep": "/opt/pirs-2.0.1/Profiles/GC-depth_Profiles/humNew.gcdep_100.dat",
    }
    expected_settings = ["PE100", "indels", "gcdep"]

    def run(self):
        # last minute settings
        self.before_run()
        try:
            rv = self.db_api.get_fastas()
            for fa in ['fasta0', 'fasta1']:
                self.module_settings[fa] = rv[fa]
        except Exception as e:
            logger.error('Pirs is missing a required modified Fasta')
            logger.error("Exception: {}".format(e))
            sys.exit(1)

        # run
        logger.info('Pirs: simulating reads ...')
        log = os.path.join(self.module_settings['outdir'], "pirs.log")
        self.module_settings['pirs_log'] = log

        cmd = "pirs simulate -l 100 -x 30 -o {outdir}" + \
            " --insert-len-mean=180 --insert-len-sd=18 --diploid " + \
            " --base-calling-profile={PE100}" + \
            " --indel-error-profile={indels}" + \
            " --gc-bias-profile={gcdep}" + \
            " --phred-offset=33 --no-gc-bias -c gzip " + \
            " -t 48 {fasta0} {fasta1} " + \
            " --no-indel-errors &> {pirs_log}"

        cmd = cmd.format(**self.module_settings)

        logger.info("pirs cmd: {}".format(cmd))
        try:
            subprocess.check_output(cmd, shell=True)
        except Exception() as e:
            logging.error('Error message %s' % e)
            raise

        logger.info('find output fastqs')
        fq1_list = glob.glob(os.path.join(self.module_settings["outdir"], 'pirs', '*1.fq.gz'))
        fq2_list = glob.glob(os.path.join(self.module_settings["outdir"], 'pirs', '*2.fq.gz'))
        self.db_api.post_reads(fq1_list[0], fq2_list[0])


###########################################################
class RSVSim(ModuleBase):
    default_settings = {}
    expected_settings = []

    def run(self):
        self.before_run()
        pass


###########################################################
class VCF2Fasta(ModuleBase):
    default_settings = {}
    expected_settings = []

    def run(self):
        self.before_run()
        pass


###########################################################
class CompositeDataset(ModuleBase):
    default_settings = {}
    expected_settings = []

    def run(self):
        self.before_run()
        pass