"""
A module defines a unit of work ( creating a truth VCF/ fasta etc )
"""

import os, sys, logging, copy
import shutil
from lib.common import PipelineExc, parse_ht_config
from lib.common import MPdb
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
        self.module_settings = {}
        self.check_and_update_input_settings()
        self.db_api = db_api
        self.lock = self.pipeline_settings["lock"]


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

    @abstractproperty
    def promised_outputs(self):
        """
        this module will create these outputs
        so next module can use this to validate expected settings
        :rtype: list
        """
        pass

    def check_and_update_input_settings(self):
        """
        for each module use the default settings
        if a settings is available in pipeline settings then override
        assert if setting is never specified
        """

        # this module promised to create these outputs
        for output in self.promised_outputs:
            self.module_settings[output] = "_PROMISED"

        for key in self.default_settings:
            self.module_settings[key] = self.default_settings[key]

        for key in self.pipeline_settings:
            # if provided use it 
            self.module_settings[key] = self.pipeline_settings[key]

        for key in self.expected_settings:
            # if neither default nor provided, then assert
            if key not in self.module_settings:
                raise PipelineExc("Module: {}, missing setting: {}".format(self.name, key))

        # workdir/ outdir module name
        for d in ["workdir", "outdir"]:
            self.module_settings[d] = os.path.join(
                self.pipeline_settings[d], self.name)

        # this step is required for init/ validation
        # we expect that the module will update this key
        for output in self.promised_outputs:
            if not self.pipeline_settings.get(output):
                self.pipeline_settings[output] = "_PROMISED"

    def create_workdir_outdir(self):
        """
        create sub directories for each module
        """
        for _d in ['workdir', 'outdir']:
            _dir = self.module_settings[_d]

            # start with clean dir
            if os.path.isdir(_dir):
                shutil.rmtree(_dir)

            try:
                self.logger.info("Create {}: {}".format(_d, _dir))
                os.makedirs(_dir)
            except Exception as e:
                msg = "Failed to create {}: {}, exception: {}".format(_d, _dir, e)
                raise PipelineExc(msg)

    def print_module_settings(self):
        self.logger.info(self.module_settings)

    def before_run(self):
        self.create_workdir_outdir()
        self.logger.info("")
        self.logger.info("- module: {}".format(self.name))
        self.exit_status = "FAILED"
        # sync module with newest pileline settings
        self.check_and_update_input_settings()

    def after_run(self):
        self.logger.info("- done: {}".format(self.name))
        self.check_and_update_promised_outputs()
        self.exit_status = "COMPLETED"
        self.update_db()
        # delete possible big variables 
        del self.module_settings


    def update_db(self):
        if hasattr(self, "update_db_keys"):
            for i in self.update_db_keys:
                try:
                    self.db_api.upload_to_db(i, self.module_settings[i])
                except:
                    raise PipelineExc("Invalid key to update db: {}".format(i))


    def check_and_update_promised_outputs(self):
        """
        make sure each module delivered on the outputs they promised
        also copy these promised outputs from the module settings to the pipeline_settings
        """

        for output in self.promised_outputs:
            if self.module_settings[output] == "_PROMISED":
                raise PipelineExc(" - missing promised result: {}".format(output))
            else:
                self.logger.info("Updating pipeline setting: {}".format(output))
                self.pipeline_settings[output] = self.module_settings[output]


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
