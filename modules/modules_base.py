"""
A module defines a unit of work ( creating a truth VCF/ fasta etc )
"""

import os
import sys
import logging
from lib.common import PipelineExc, parse_ht_config
from abc import ABCMeta, abstractmethod, abstractproperty

logger = logging.getLogger(__name__)
this_dir_path = os.path.dirname(os.path.realpath(__file__))


###########################################################
class ModuleBase(object):
    """
    Base class for all modules (pipelines)
    """
    __metaclass__ = ABCMeta

    def __init__(self, pipeline_settings, db_api):
        self.name = self.__class__.__name__
        logger.info("- validating module: {}".format(self.name))
        self.pipeline_settings = pipeline_settings
        self.dataset_name = None
        self.module_settings = None
        self.parse_settings()
        self.db_api = db_api

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
        self.module_settings['outdir'] = os.path.join(
            self.outdir, self.pipeline_settings['dataset_name'], self.name)
        if not os.path.isdir(self.module_settings['outdir']):
            try:
                logger.info("create module outdir: {}".format(self.module_settings['outdir']))
                os.makedirs(self.module_settings['outdir'])
            except Exception as e:
                msg = "Failed to create directory: {}, exception: {}".format(
                    self.module_settings['outdir'], e)
                raise PipelineExc(msg)

        self.validate_module_settings()

    def print_module_settings(self):
        logger.info(self.module_settings)

    def before_run(self):
        logger.info("")
        logger.info("- module: {}".format(self.name))

    def after_run(self):
        logger.info("- done: {}".format(self.name))

    def add_module_settings_to_saved_outputs_dict(self):
        for key in self.module_settings:
            self.db_api.outputs_dict[key] = self.module_settings[key]

    @abstractmethod
    def run(self):
        pass

    def get_dataset_ref(self):
        rv = self.db_api.get_dataset_ref_info()

        keys = ["ref_type", "fasta_file"]
        for key in keys:
            self.module_settings[key] = rv[key]
            if not self.module_settings[key]:
                msg = "{}, dataset {} missing HT reference file: {}"
                msg = msg.format(self.name, self.dataset_name, key)
                raise PipelineExc(msg)

    def get_ht_cfg(self):
        key = "hash_table5"
        hash_dir = self.get_from_db(self.dataset_name, key)
        if hash_dir[:-1] = "\\":
            hash_dir = hash_dir[:-1]
        cfg_file = os.path.join(hash_dir, "hash_table.cfg")
        try:
            res = parse_ht_config(cfg_file)
            return res
        except:
            raise PipelineExc("Failed to parse the HT config file for contig lengths")

    def get_target_bed(self):
        self.module_settings["target_bed"] = self.db_api.get_target_bed()

    def post_fastas(self, fasta0, fasta1):
        self.db_api.set_fastas(fasta0, fasta1)

    @staticmethod
    def fail(msg, mod):
        msg = "Invalid settings for module: {}, {}".format(mod, msg)
        logger.error(msg)
        sys.exit(1)
