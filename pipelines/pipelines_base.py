from abc import ABCMeta, abstractproperty
from lib.db_api import DBAPI
from lib.common import PipelineExc, MPdb
import logging, sys, os, copy

logger = logging.getLogger(__name__)


###########################################################
# all pipelines must inherit from this base

class PipelinesBase(object):
    """
    Base class for all other classes in this module
    """
    __metaclass__ = ABCMeta

    def __init__(self, settings, pipeline_idx, simulator_settings):
        self._class_name = self.__class__.__name__
        self.simulator_settings = simulator_settings
        self.pipeline_idx = pipeline_idx
        self.validate_pipeline_settings(settings)
        self.update_settings_for_multiple_runs()
        self.create_work_and_outdir()
        self.dataset_name = self.pipeline_settings["dataset_name"]
        self.db_api = DBAPI(self.dataset_name)
        self.module_instances = []
        self.db_api.dataset_create_or_update(
            self.pipeline_settings["dataset_name"],
            self.pipeline_settings["reference"])


    def validate_pipeline_settings(self, settings):
        """
        validate pipeline settings
        """
        # before we have a local logger, we log globally 
        msg = "Validating pipeline settings: {}".format(self.__class__.__name__)
        logger.info(msg)

        # avoid sharing settings with other runs of this pipeline
        self.pipeline_settings = copy.deepcopy(settings)

        # required settings for each pipeline
        expected = ["outdir", "workdir", "dataset_name", "reference"] 

        for i in expected:
            if not i in self.pipeline_settings:
                self.logger.error("\nMissing pipeline setting: {}".format(i))
                sys.exit(1)

    def update_settings_for_multiple_runs(self):
        """ 
        create subdirs for each dataset
        """
        if self.pipeline_idx > 0:
            self.pipeline_settings["dataset_name"] += "_{}".format(
                self.pipeline_idx) 

        for d in ["workdir", "outdir"]:
            self.pipeline_settings[d] = os.path.join(
                self.pipeline_settings[d],
                self.pipeline_settings["dataset_name"])


    def create_work_and_outdir(self):
        """
        need to create the outdir before we can start using the local logger
        """
        for _d in ['workdir', 'outdir']:
            _dir = self.pipeline_settings[_d]

            if not os.path.isdir(_dir):
                try:
                    os.makedirs(_dir)
                except Exception as e:
                    msg = "Failed to create {}: {}, exception: {}".format(
                        _d, _dir, e)
                    raise PipelineExc(msg)        

    @property
    def name(self):
        return self._class_name

    @abstractproperty
    def modules(self):
        """ :rtype: list """
        pass

    def instantiate_modules(self):
        """ validate settings prior to run """
        for C in self.modules:
            inst = C(self.pipeline_settings, self.db_api, self.logger, 
                     self.simulator_settings)
            self.module_instances.append(inst)

    def run(self):
        self.logger.info("\nPIPELINE: {}".format(self.name))
        # are there any values we just need to upload to the dataset db?
        # these keys are not intended to be used be modules
        for_upload = self.pipeline_settings.get("for_upload", None)
        if for_upload:
            for key in for_upload:
                self.db_api.upload_to_db(key, for_upload[key])
                
        # run modules
        for inst in self.module_instances:
            inst.before_run()
            inst.run()
            inst.after_run()




        
