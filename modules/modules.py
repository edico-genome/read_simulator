"""
A module defines a unit of work ( creating a truth VCF/ fasta etc )
"""

import os
import sys
import glob
import logging
import vlrd_lib
import subprocess
from lib.db_api import DBAPI
from abc import ABCMeta, abstractmethod, abstractproperty

logger = logging.getLogger(__name__)


###########################################################
class ModuleBase(object):
    """
    Base class for all modules (pipelines)
    """
    __metaclass__ = ABCMeta

    @abstractproperty
    def name(self):
        pass

    @abstractproperty
    def default_settings(self):
        """
        a dictionary of key: value pairs that will be
        used as default module settings
        :rtype: dict
        """
        pass

    @abstractproperty
    def expected_settings(self):
        """
        list the keys that should be in module settings
        used to validate module settings prior to run
        :rtype: list
        """
        pass

    def __init__(self, pipeline_settings):
        self.pipeline_settings = pipeline_settings
        self.dataset_name = None
        self.module_settings = None
        self.parse_settings()
        self.db_api = DBAPI(self.dataset_name)

    def validate_module_settings(self):
        for key in self.default_settings:
            if key not in self.module_settings:
                self.module_settings[key] = self.default_settings[key]
        for key in self.expected_settings:
            assert self.module_settings[key], "missing key {} for {}".format(key, self.name)

    def parse_settings(self):
        for key in ["dataset_name", "module_settings", "outdir", "reference"]:
            value = self.pipeline_settings.get(key, None)
            if not value:
                self.fail("'{}' not specified".format(key), self.name)
            self.__setattr__(key, value)

        self.module_settings = self.pipeline_settings["module_settings"].get(self.name, {})
        if not self.module_settings:
            logging.warning("'{}' not specified within run_config: module_settings".format(self.name))

        # create sub directory for each module
        self.module_settings['outdir'] = os.path.join(self.outdir, self.name)
        if not os.path.isdir(self.module_settings['outdir']):
            try:
                os.makedirs(self.module_settings['outdir'])
            except Exception as e:
                logger.error("Failed to create directory: {}, exception: {}".
                             format(self.module_settings['outdir'], e))

        self.validate_module_settings()

    def print_module_settings(self):
        logger.info(self.module_settings)

    @abstractmethod
    def run(self):
        pass

    def get_dataset_ref(self):
        rv = self.db_api.get_dataset_ref_info()
        self.module_settings["ref_type"] = rv["ref_type"]
        self.module_settings["ref_fasta"] = rv["fasta"]

    def get_target_bed(self):
        self.module_settings["target_bed"] = self.db_api.get_target_bed()

    def post_fastas(self, fasta0, fasta1):
        self.db_api.set_fastas(fasta0, fasta1)

    @staticmethod
    def fail(msg, mod):
        msg = "Invalid settings for module: {}, {}".format(mod, msg)
        logger.error(msg)
        sys.exit(1)


###########################################################
class VLRDmod(ModuleBase):
    name = "vlrd"
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
        res = vlrd_lib.create_truth_vcf_and_fastas(self.module_settings)
        logger.info(res)
        self.db_api.set_fastas(res['fasta0'], res['fasta1'])
        self.db_api.post_truth_vcf(res['truth_vcf'])


###########################################################
class Pirs(ModuleBase):
    name = "pirs"
    default_settings = {
        "PE100": "/opt/pirs-2.0.1/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz",
        "indels": "/opt/pirs-2.0.1/Profiles/InDel_Profiles/phixv2.InDel.matrix",
        "gcdep": "/opt/pirs-2.0.1/Profiles/GC-depth_Profiles/humNew.gcdep_100.dat",
    }
    expected_settings = ["PE100", "indels", "gcdep"]

    def run(self):
        # last minute settings
        rv = self.db_api.get_fastas()
        self.module_settings['fasta0'] = rv['fasta0']
        self.module_settings['fasta1'] = rv['fasta1']
        for f in ['fasta0', 'fasta1']:
            assert os.path.isfile(self.module_settings[f]), \
                'Modified %s is not a valid file' % f

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
    name = "rsvsim"
    default_settings = {}
    expected_settings = []

    def run(self):
        pass


###########################################################
class VCF2Fasta(ModuleBase):
    name = 'vcf2fasta'
    default_settings = {}
    expected_settings = []

    def run(self):
        logger.info('VCF2Fasta: generating truth fasta ...')


###########################################################
class CompositeDataset(ModuleBase):
    name = "composite_dataset"
    default_settings = {}
    expected_settings = []

    def run(self):
        logger.info('CompositeDataset: combining datasets with specified allele frequencies ...')
