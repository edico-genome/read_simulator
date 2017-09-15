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
from lib.common import PipelineExc

logger = logging.getLogger(__name__)
this_dir_path = os.path.dirname(os.path.realpath(__file__))


###########################################################
def check_output(cmd):
    """ helper function to run cmd """
    logger.debug("Run cmd: {}".format(cmd))
    res = subprocess.check_output(cmd, shell=True)
    logger.debug("Response: {}".format(res))
    return res


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

        # use case 1: is one contigs primary & the other contig an alt?
        elif contig_is_primary["contig-1"] != contig_is_primary["contig-2"]:
            use_case = 1
            # make contig-2 the primary
            if contig_is_primary["contig-1"]:
                _tmp = self.module_settings["contig-1"]
                self.module_settings["contig-1"] = self.module_settings["contig-2"]
                self.module_settings["contig-2"] = _tmp

        # use case 1/ 2: are both alts?
        else:
            # are they the same alt?
            if self.module_settings["contig-1"] == self.module_settings["contig-2"]:
                use_case = 2
            else:
                use_case = 3

        self.module_settings["contigs_combo_type"] = use_case
        assert self.module_settings["contigs_combo_type"] in [1, 2, 3],\
            "invalid contigs specified"

    def get_contig_ranges(self):
        """ lookup alt contig length in dict """
        logger.info("Get contig ranges")

        with open(self.module_settings["dict_file"], 'r') as stream:
            # set from equal to 1
            self.module_settings["contig-1-from"] = 1
            self.module_settings["contig-2-from"] = 1

            # set "to" equal to length of contig
            for line in stream:
                len_col = line.split()[2]
                if "LN" in len_col:
                    for c_idx in ["contig-1", "contig-2"]:
                        c_idx_to = c_idx + "-to"
                        contig = self.module_settings[c_idx]
                        if contig in line:
                            try:
                                c_len = int(len_col.replace("LN:", ""))
                                self.module_settings[c_idx_to] = c_len
                                logger.debug("{}: {}".format(c_idx, c_len))
                            except Exception as e:
                                logger.error("Contig len not found in line: {}".format(line))
                                logger.error(e)

        for i in ["contig-1-from", "contig-2-from", "contig-1-to", "contig-2-to"]:
            assert isinstance(self.module_settings[i], int), \
                "Invalid contig ranges detected"

    def get_subrange_of_primary(self):
        """ find _from and _to indexes in primary that corresponds to the contig 1
         contig 1 will always be an alt contig """
        logger.info("find indexes in primary that corresponds to alt-contig 1 ")
        _alt_to_pri_script = os.path.join(this_dir_path, "alt_contig", "convert_alt_to_pri_coords.pl")
        cmd = "{} {} {} 1 {} 2> /dev/null".format(
            _alt_to_pri_script,
            self.module_settings["alt-sam"],
            self.module_settings["contig-1"],
            self.module_settings["contig-1-to"])

        cmd = cmd.format(**self.module_settings)
        res = check_output(cmd)

        try:
            _contig, _from, _to = res.split()[0:3]
            self.module_settings["contig-primary-from"] = _from
            self.module_settings["contig-primary-to"] = _to
        except Exception as e:
            logger.error('Failed to align alt contig 1 to primary')
            logger.error(e)
            sys.exit(1)

    def get_subrange_of_contig2(self):
        """ find _from and _to indexes in alt-contig 2 that corresponds to subrange in primary """
        logger.info("find indexes in alt-contig 2 that corresponds to subrange in primary ")
        _pri_to_alt_script = os.path.join(this_dir_path, "alt_contig", "convert_pri_to_alt_coords.pl")

        cmd = "{} {} {} {} {} 2> /dev/null".format(
            _pri_to_alt_script,
            self.module_settings["alt-sam"],
            self.module_settings["contig-2"],
            self.module_settings["contig-primary-from"],
            self.module_settings["contig-primary-to"])

        res = check_output(cmd)
        try:
            _contig, _from, _to = res.split()[0:3]
            self.module_settings["contig-2-from"] = _from
            self.module_settings["contig-2-to"] = _to
        except Exception as e:
            logger.error('Failed to align alt contig 1 to primary')
            logger.error(e)
            sys.exit(1)

    def get_ranges_where_contigs_match(self):
        """ get coordinates which we'll use for the truth vcf generation
            the coordinates depend on the combo types
            1. primary + alt_contig ( we also want to standardize this order )
            2. alt_contig1 + alt_contig1
            3. alt_contig1 + alt_contig2(any 2 different, but related, contigs)
        """

        # contig-1 will always be an alt contig ( use case 1, 2 and 3 )
        self.get_subrange_of_primary()

        # if contig2 is also an alt contig, and differs from contig 1
        # then we need to update it's range so that it aligns to contig 1
        if self.module_settings["contigs_combo_type"] == 3:
            self.get_subrange_of_contig2()

    def create_modified_fastas(self):
        fastas = {"contig-1": "{outdir}/mod_fasta_0.fa".format(**self.module_settings),
                  "contig-2": "{outdir}/mod_fasta_1.fa".format(**self.module_settings)}

        for c in ["contig-1", "contig-2"]:
            cmd = "samtools faidx {} {}:{}-{} > {}"
            cmd = cmd.format(self.module_settings["fasta_file"],
                             self.module_settings[c],
                             self.module_settings[c + "-from"],
                             self.module_settings[c + "-to"],
                             fastas[c])
            try:
                subprocess.check_call(cmd, shell=True)
            except Exception as e:
                raise PipelineExc(e)

        self.db_api.set_fastas(fastas["contig-1"], fastas["contig-2"])

    def create_truth_vcf(self):
        script_path = os.path.join(this_dir_path, "alt_contig", "fasta_sam_to_vcf.pl")
        self.module_settings["truth_vcf"] = os.path.join(self.module_settings["outdir"], "alt_contig_truth.vcf")
        cmd = script_path + " {fasta_file} {alt-sam} {contig-1}:{contig-1-from}-{contig-1-to} " + \
            "{contig-2}:{contig-2-from}-{contig-2-to} > {truth_vcf}"
        cmd = cmd.format(**self.module_settings)
        logger.info("create truth vcf cmd: {}".format(cmd))

        try:
            check_output(cmd)
        except Exception as e:
            raise Exception('Failed to create truth VCF: {}'.format(e))
        self.db_api.post_truth_vcf(self.module_settings["truth_vcf"])

    def run(self):
        self.get_dataset_ref()
        self.identify_contig_types_and_use_case()
        self.get_contig_ranges()
        self.get_ranges_where_contigs_match()
        self.create_modified_fastas()
        self.create_truth_vcf()
        # self.add_module_settings_to_saved_outputs_dict()


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

    def update_run_settings(self):
        """ get settings from previous runs """
        for key in ["contig", "alt-contig-1", "alt-contig-2", "alt-sam", "read_info"]:
            if key in self.db_api.outputs_dict:
                self.module_settings[key] = self.db_api.outputs_dict[key]

    def run(self):
        self.module_settings["truth_sam"] = \
            os.path.join(self.module_settings["outdir"], "truth.sam")
        self.update_run_settings()

        # test for homozygous alts
        if self.module_settings["alt-contig-1"] == self.module_settings["alt-contig-2"]:
            script_path = os.path.join(this_dir_path, "alt_contig", "xform_pirs_read_info.pl")
            self.module_settings["xform_pirs_read_info"] = script_path
            cmd = "{xform_pirs_read_info} <(zcat {read_info}) {alt-sam} > {truth_sam}".format(**self.module_settings)
        else:
            raise PipelineExc("Use case not supported")

        try:
            res = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
            logging.info("Response: {}".format(res))
            logging.info("Truth sam created: {}".
                         format(self.module_settings["truth_sam"]))
            self.db_api.upload_to_db('sam_gold', self.module_settings["truth_sam"])
        except Exception as e:
            logging.error('Error message %s' % e)


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
