"""
A module defines a unit of work ( creating a truth VCF/ fasta etc )
"""

import os
import sys
import glob
import logging
import subprocess
from modules_base import ModuleBase
from lib.common import PipelineExc
from enum import Enum
from vlrd import vlrd_functions


logger = logging.getLogger(__name__)
this_dir_path = os.path.dirname(os.path.realpath(__file__))


###########################################################
def check_output(cmd):
    """ helper function to run cmd """
    logger.info("Run cmd: {}".format(cmd))
    res = subprocess.check_output(cmd, shell=True)
    logger.info("Response: {}".format(res))
    return res


###########################################################
class CNVExomeModFastas(ModuleBase):
    default_settings = {}
    expected_settings = [
        "nHomozygousDeletions",
        "nHeterozygousDeletions",
        "nTandemDuplications", "maxEventLength",
        "minEventLength", "outdir", "fa_file",
        "cnv_db", "target_chrs", "target_bed"]

    def get_refs(self):
        self.get_dataset_ref()
        self.get_target_bed()
        for key in ["ref_type", "fasta_file"]:
            if not self.module_settings[key]:
                msg = "{}, dataset {} missing: {}".format(self.name, self.dataset_name, key)
                raise PipelineExc(msg)

    def run(self):
        self.get_refs()
        _mod_fastas_script = os.path.join(this_dir_path, "cnv_exomes", "RSVSim_generate_modified_genome.R")

        cmd =_mod_fastas_script
        cmd += " --nHomozygousDeletions {nHomozygousDeletions}"
        cmd += " --nHeterozygousDeletions {nHeterozygousDeletions}"
        cmd += " --nTandemDuplications {nTandemDuplications}"
        cmd += " --maxEventLength {maxEventLength}"
        cmd += " --minEventLength {minEventLength}"
        cmd += " --outdir {outdir}"
        cmd += " --fa_file {fasta_file}"
        cmd += " --cnv_db {cnv_db}"
        cmd += " --target_chrs {target_chrs}"
        cmd += " --target_bed {target_bed}"
        cmd = cmd.format(**self.module_settings)

        try:
            logger.info(cmd)
            subprocess.check_call(cmd, shell=True)
        except Exception as e:
            raise PipelineExc(e)

        fastas = {"contig-1": "{outdir}/genome_rearranged1.fasta".format(**self.module_settings),
                  "contig-2": "{outdir}/genome_rearranged2.fasta".format(**self.module_settings)}

        csv_files = ["insertions1.csv",
                     "insertions2.csv",
                     "deletions.csv",
                     "deletions1.csv",
                     "deletions2.csv",
                     "tandemDuplications.csv",
                     "tandemDuplications1.csv",    
                     "tandemDuplications2.csv"]

        self.db_api.set_fastas(fastas["contig-1"], fastas["contig-2"])


###########################################################
class VLRDVCF(ModuleBase):
    default_settings = {}
    expected_settings = ["varrate"]

    def get_refs(self):
        self.get_dataset_ref()
        self.get_target_bed()
        for key in ["ref_type", "fasta_file", "dict_file"]:
            if not self.module_settings[key]:
                logger.info("{}, dataset {} missing: {}".format(
                    self.name, self.dataset_name, key))
                sys.exit(1)

    def run(self):
        self.get_refs()
        res = vlrd_functions.create_truth_vcf_and_fastas(self.module_settings)
        logger.info(res)
        self.db_api.set_fastas(res['fasta0'], res['fasta1'])
        self.db_api.post_truth_vcf(res['truth_vcf'])


###########################################################
class ContigComboType(Enum):
    alt_prim = 1
    same_alts = 2
    diff_alts = 3


class AltContigVCF(ModuleBase):
    default_settings = {}
    expected_settings = ["contig-1", "contig-2", "alt-sam"]

    def identify_contig_types_and_use_case(self):
        """
        for the 2 contigs; the user can specify any of the following 3 contig combo cases
        1. primary + alt_contig ( we also want to standardize this order )
        2. alt_contig1 + alt_contig1
        3. alt_contig1 + alt_contig2(any 2 different, but related, contigs)
        """

        # identify type for each contig
        contig_is_primary = {}
        for contig in ["contig-1", "contig-2"]:
            contig_is_primary[contig] = True
            for sufx in ['_alt', 'HLA']:  # these markers help identify alt contigs
                if sufx in self.module_settings[contig]:
                    contig_is_primary[contig] = False

        # both contigs should not be primary
        if all(contig_is_primary.values()):  # if all([contig_is_primary[key] for key in ["contig-1", "contig-2"]):
            raise PipelineExc("Use case not supported where both contigs are primary")

        # one contig is an alt and the other is a primary
        elif contig_is_primary["contig-1"] != contig_is_primary["contig-2"]:
            self.module_settings["contigs_combo_type"] = ContigComboType.alt_prim

            # standardize with contig-2 being the primary
            # contig 1 will now always be an alt
            if contig_is_primary["contig-1"]:
                _tmp = self.module_settings["contig-1"]
                self.module_settings["contig-1"] = self.module_settings["contig-2"]
                self.module_settings["contig-2"] = _tmp

        # both contigs are alt-contigs
        else:
            # are they the same alt?
            if self.module_settings["contig-1"] == self.module_settings["contig-2"]:
                self.module_settings["contigs_combo_type"] = ContigComboType.same_alts
            else:
                self.module_settings["contigs_combo_type"] = ContigComboType.diff_alts

    def get_contig_length_from_cfg(self, contig_key):
        """ lookup alt contig length in dict """
        assert contig_key in ["contig-1", "contig-2"]
        logger.info("Get contig {} ranges".format(contig_key))

        try:
            c_idx = contig_key
            c_idx_from = c_idx + "-from"
            c_idx_to = c_idx + "-to"
            contig_name = self.module_settings[c_idx]
            self.module_settings[c_idx_from] = 1
            self.module_settings[c_idx_to] = self.module_settings["ht_table_config"][contig_name]
        except Exception as e:
            raise PipelineExc("Contig not found: {}".format(e))

        for i in [c_idx_from, c_idx_to]:
            assert isinstance(self.module_settings[i], int), \
                "Invalid contig ranges detected"

    def get_contig_start_stop_indexes(self):
        """ a specific range will be available on the primary, alt 1 and alt 2 
        if assume the whole of alt 1 will be represented
        now find corresponding region in primary """

        # get contig 1 start stop indexes
        self.get_contig_length_from_cfg("contig-1")

        # get primary contig start stop indexes
        logger.info("find primary contig start stop indexes")
        _alt_to_pri_script = os.path.join(this_dir_path, "alt_contig", "convert_alt_to_pri_coords.pl")
        cmd = "{} {} {} 1 {} 2> /dev/null".format(
            _alt_to_pri_script,
            self.module_settings["alt-sam"],
            self.module_settings["contig-1"],
            self.module_settings["contig-1-to"])

        cmd = cmd.format(**self.module_settings)
        res = check_output(cmd)

        try:
            _from, _to = res.split()[1:3]
            self.module_settings["contig-primary-from"] = _from
            self.module_settings["contig-primary-to"] = _to
        except Exception as e:
            raise PipelineExc('Failed to find start-stop index in primary {}'.format(e))

        # get contig 2 start stop indexes ( context specific )
        logger.info("find indexes in alt-contig 2 that corresponds to subrange in primary ")
        if self.module_settings["contigs_combo_type"] == ContigComboType.alt_prim:
            self.module_settings["contig-2-from"] = self.module_settings["contig-primary-from"]
            self.module_settings["contig-2-to"] = self.module_settings["contig-primary-to"]
        elif self.module_settings["contigs_combo_type"] == ContigComboType.same_alts:
            self.module_settings["contig-2-from"] = self.module_settings["contig-1-from"] 
            self.module_settings["contig-2-to"] = self.module_settings["contig-1-to"] 
        elif self.module_settings["contigs_combo_type"] == ContigComboType.diff_alts:
            # find range in contig 2 ( alt contig ) that maps to alt contig 1
            _pri_to_alt_script = os.path.join(this_dir_path, "alt_contig", "convert_pri_to_alt_coords.pl")
            cmd = "{} {} {} {} {} 2> /dev/null".format(
                _pri_to_alt_script,
                self.module_settings["alt-sam"],
                self.module_settings["contig-2"],
                self.module_settings["contig-primary-from"],
                self.module_settings["contig-primary-to"])
            res = check_output(cmd)
            try:
                _from, _to = res.split()[4:6]
                self.module_settings["contig-2-from"] = _from
                self.module_settings["contig-2-to"] = _to
            except Exception as e:
                raise PipelineExc('Failed to find start-stop indexes in alt-contig-2 {}'.format(e))

    def create_modified_fastas(self):
        """ create fastas for read simulator """

        def check_call():
            try:
                logger.info(cmd)
                subprocess.check_call(cmd, shell=True)
            except Exception as e:
                raise PipelineExc(e)

        fastas = {"contig-1": "{outdir}/mod_fasta_0.fa".format(**self.module_settings),
                  "contig-2": "{outdir}/mod_fasta_1.fa".format(**self.module_settings)}

        for c in ["contig-1", "contig-2"]:
            cmd = "samtools faidx {} {}:{}-{} > {}"
            cmd = cmd.format(self.module_settings["fasta_file"],
                             self.module_settings[c],
                             self.module_settings[c + "-from"],
                             self.module_settings[c + "-to"],
                             fastas[c])
            check_call()

        self.db_api.set_fastas(fastas["contig-1"], fastas["contig-2"])

    def create_truth_vcf(self):
        script_path = os.path.join(this_dir_path, "alt_contig", "fasta_sam_to_vcf.pl")
        self.module_settings["truth_vcf"] = os.path.join(self.module_settings["outdir"], "alt_contig_truth.vcf")
        cmd = script_path + " {fasta_file} {alt-sam} {contig-1}:{contig-1-from}-{contig-1-to} " + \
            "{contig-2}:{contig-2-from}-{contig-2-to} > {truth_vcf}"
        cmd = cmd.format(**self.module_settings)
        try: 
            check_output(cmd)
        except Exception as e:
            raise PipelineExc("Failed to create truth VCF: {}".format(e))
        self.db_api.post_truth_vcf(self.module_settings["truth_vcf"])

    def run(self):
        self.get_dataset_ref()
        self.get_ht_cfg()
        self.identify_contig_types_and_use_case()
        self.get_contig_start_stop_indexes()
        self.create_modified_fastas()
        self.create_truth_vcf()
        self.add_module_settings_to_saved_outputs_dict()


###########################################################
class Capsim(ModuleBase):
    default_settings = {}
    expected_settings = ["number-of-reads", "read-len", "fragment-size"]

    def run(self):
        # get modified fastas
        try:
            rv = self.db_api.get_fastas()
            for f in ['fasta0', 'fasta1']:
                self.module_settings[f] = rv[f]
        except Exception as e:
            raise PipelineExc('Capsim is missing a required modified Fasta: {}'.format(e))

        # run
        logger.info('Capsim: simulating reads ...')
        script_path = os.path.join(this_dir_path, "modules", "cnv_exomes", "run_capsim.sh")
        log = os.path.join(self.module_settings['outdir'], "capsim.log")
        self.module_settings['capsim_log'] = log

        cmd = script_path + \
            " -n {number-of-reads}" + \
            " -l {read-len}" + \
            " -f {fragment-size}" + \
            " -o {outdir}" + \
            " -1 {fasta0} -2 {fasta1}"

        cmd = cmd.format(**self.module_settings)

        logger.info("Capsim cmd: {}".format(cmd))
        try:
            subprocess.check_output(cmd, shell=True)
        except Exception() as e:
            logging.error('Error message %s' % e)
            raise

        # update results
        logger.info('find the Capsim generated FQs and upload to DB')
        fq_lists = []
        for i in ["1", "2"]:
            p = os.path.join(self.module_settings["outdir"], 'output_' + i + '.fastq.gz')
            fq_lists.append(glob.glob(p))
        self.db_api.post_reads(fq_lists[0][0], fq_lists[1][0])


###########################################################
class Pirs(ModuleBase):
    default_settings = {
        "PE100": "/opt/pirs-2.0.1/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz",
        "indels": "/opt/pirs-2.0.1/Profiles/InDel_Profiles/phixv2.InDel.matrix",
        "gcdep": "/opt/pirs-2.0.1/Profiles/GC-depth_Profiles/humNew.gcdep_100.dat",
    }
    expected_settings = ["PE100", "indels", "gcdep"]

    def run(self):
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

        cmd = "pirs simulate -l 100 -x 30 -o {outdir}/pirs" + \
              " --insert-len-mean=180 --insert-len-sd=18 --diploid " + \
              " --base-calling-profile={PE100}" + \
              " --indel-error-profile={indels}" + \
              " --gc-bias-profile={gcdep}" + \
              " --phred-offset=33 --no-gc-bias -c gzip " + \
              " -t 48 {fasta0} {fasta1} " + \
              " --no-indel-errors &> {pirs_log}"

        cmd = cmd.format(**self.module_settings)

        logger.info("Pirs cmd: {}".format(cmd))
        try:
            subprocess.check_output(cmd, shell=True)
        except Exception() as e:
            logging.error('Error message %s' % e)
            raise

        logger.info('find the pirs generated FQs and upload to DB')
        fq_lists = []
        for i in ["1", "2"]:
            p = os.path.join(self.module_settings["outdir"], 'pirs*' + i + '.fq.gz')
            fq_lists.append(glob.glob(p))
        self.db_api.post_reads(fq_lists[0][0], fq_lists[1][0])

        p = os.path.join(self.module_settings["outdir"], 'pirs*read.info.gz')
        l = glob.glob(p)
        if os.path.isfile(l[0]):
            self.db_api.outputs_dict['read_info'] = l[0]
        else:
            raise PipelineExc("Failed to find read_info: {}".format(p))


###########################################################
class AltContigPirsTruthSam(ModuleBase):
    default_settings = {}
    expected_settings = []

    def run(self):
        # truth sam
        self.module_settings["truth_sam"] = \
            os.path.join(self.module_settings["outdir"], "truth.sam")

        # xform path
        xform_script_path = os.path.join(this_dir_path, "alt_contig", "xform_pirs_read_info.pl")

        # get settings from previous runs
        for key in ["contig", "contig-1", "contig-2", "alt-sam", "read_info",
                    "contigs_combo_type", "contig-1-from", "contig-1-to", "contig-2-from",
                    "contig-2-to"]:
            if key in self.db_api.outputs_dict:
                self.module_settings[key] = self.db_api.outputs_dict[key]

        bcf_modified_name = {}
        for i in ["1", "2"]:
            _contig_name = self.module_settings["contig-{}".format(i)],
            _from = self.module_settings["contig-{}-from".format(i)],
            _to = self.module_settings["contig-{}-to".format(i)]

            bcf_modified_name[_contig_name] = "{}:{}-{}".format(_contig_name, _from, _to)

        options = {
            "xform": xform_script_path,
            "read_info": self.module_settings["read_info"],
            "alt-sam": self.module_settings["alt-sam"],
            "bcf-contig-1": bcf_modified_name["contig-1"],
            "bcf-contig-2": bcf_modified_name["contig-2"],
            "orig-contig-1": self.module_settings["contig-1"],
            "orig-contig-2": self.module_settings["contig-2"],
            "truth_sam": self.module_settings["truth_sam"]
        }

        cmd = "{xform} --read-info <(zcat {read_info}) --alt-sam {alt-sam}"
        cmd += " --contig1 {bcf-contig-1} --contig1-rename {orig-contig-1} --contig1-offset {contig-1-from}"
        cmd += " --contig2 {bcf-contig-2} --contig2-rename {orig-contig-2} --contig2-offset {contig-2-from}"
        cmd += " > {truth_sam}"

        cmd = cmd.format(**options)

        try:
            res = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
            logging.info("Response: {}".format(res))
            logging.info("Truth sam created: {}".format(self.module_settings["truth_sam"]))
            self.db_api.upload_to_db('sam_gold', self.module_settings["truth_sam"])
        except Exception as e:
            raise PipelineExc('Failed to create gold Sam. Error message: %s' % e)


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
