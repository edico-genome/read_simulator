"""
A module defines a unit of work ( creating a truth VCF/ fasta etc )
"""

import os, sys, logging, copy
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

    def __init__(self, pipeline_settings, db_api, local_logger):
        self.name = self.__class__.__name__
        self.logger = local_logger
        self.logger.info("- validating module: {}".format(self.name))
        self.pipeline_settings = pipeline_settings
        self.module_settings = None
        self.validate_module_settings()
        self.create_workdir_outdir()
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
        """
        for each module use the default settings
        if a settings is available in pipeline settings then override
        assert if setting is never specified
        """
        self.module_settings = copy.deepcopy(self.default_settings)
        for key in self.expected_settings:
            if key in self.pipeline_settings:
                self.module_settings[key] = self.pipeline_settings[key]
            if key not in self.module_settings:
                msg = "Module: {}, missing settings: {}".format(self.name, key)
                raise PipelineExc(msg)


    def create_workdir_outdir(self):
        """
        create sub directories for each module
        """
        # import pdb; pdb.set_trace()

        # workdir 
        self.module_settings['workdir'] = os.path.join(
             self.pipeline_settings['workdir'],
             self.pipeline_settings['dataset_name'], self.name)

        # outdir
        self.module_settings['outdir'] = os.path.join(
            self.pipeline_settings['outdir'],
            self.pipeline_settings['dataset_name'], self.name)

        for _d in ['workdir', 'outdir']:
            _dir = self.module_settings[_d]
            if not os.path.isdir(_dir):
                try:
                    self.logger.info("Create {}: {}".format(_d, _dir))
                    os.makedirs(_dir)
                except Exception as e:
                    msg = "Failed to create {}: {}, exception: {}".format(
                        _d, _dir, e)
                    raise PipelineExc(msg)

    def print_module_settings(self):
        self.logger.info(self.module_settings)

    def before_run(self):
        self.logger.info("")
        self.logger.info("- module: {}".format(self.name))
        self.exit_status = "FAILED"

    def after_run(self):
        self.logger.info("- done: {}".format(self.name))
        self.exit_status = "COMPLETED"

    def add_module_settings_to_saved_outputs_dict(self):
        for key in self.module_settings:
            self.db_api.outputs_dict[key] = self.module_settings[key]

    @abstractmethod
    def run(self):
        """ need to be overriden in each subclass """
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
        hash_dir = self.db_api.get_from_db(self.dataset_name, key)
        if hash_dir[:-1] == "\\":
            hash_dir = hash_dir[:-1]
        self.module_settings["hash_table_cfg"] = os.path.join(hash_dir, "hash_table.cfg")
        try:
            res = parse_ht_config(cfg_file)
            return res
        except:
            raise PipelineExc("Failed to parse the HT config file for contig lengths")

    def get_target_bed(self):
        self.module_settings["target_bed"] = self.db_api.get_target_bed()

    def post_fastas(self, fasta0, fasta1):
        self.db_api.set_fastas(fasta0, fasta1)
