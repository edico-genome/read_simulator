from abc import ABCMeta, abstractproperty
from lib.db_api import DBAPI
from lib.common import PipelineExc
import logging, sys, os

logger = logging.getLogger(__name__)


###########################################################
# all pipelines must inherit from this base

class PipelinesBase(object):
    """
    Base class for all other classes in this module
    """
    __metaclass__ = ABCMeta

    def expected_inputs(self):
        expected = ["workdir"] 
        for i in expected:
            if not i in self.pipeline_settings:
                raise PipelineExc("Missing setting: {}".format(i))

    def __init__(self, settings):
        msg = "validating pipeline: {}".format(self.__class__.__name__)
        logger.info(msg)
        self.exit_status = "FAILED"
        self.pipeline_settings = settings
        self.expected_inputs()
        self.dataset_name = self.pipeline_settings["dataset_name"]
        self.db_api = DBAPI(self.dataset_name)
        self.module_instances = []
        self.instantiate_modules()
        self._class_name = self.__class__.__name__
        self.db_api.dataset_create_or_update(
            self.pipeline_settings["dataset_name"],
            self.pipeline_settings["reference"])

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
            #  all modules in same pipeline should share a db instance
            inst = C(self.pipeline_settings, self.db_api)
            self.module_instances.append(inst)

    def run(self):

        # log file per pipeline
	log_path = os.path.join(
	    self.pipeline_settings['outdir'],
            self.pipeline_settings['dataset_name'],
            self.pipeline_settings['pipeline_name']+'.log')
        fh = logging.FileHandler(log_path)

        logger.addHandler(fh)
        logger.info("\nPIPELINE: {}".format(self.name))
        for inst in self.module_instances:
            inst.logger.addHandler(fh)
            try:     
                inst.before_run()
                inst.run()
                inst.after_run()
                self.exit_status = "COMPLETED"
            except PipelineExc as e:
                msg = "Pipeline Failed: {}; reason: {}.\n\nContinue with next pipeline ..."
                msg = msg.format(self.name, e)
                self.exit_status = "FAILED"
          	raise PipelineExc(msg)
            except Exception as e:
                msg = "Pipeline Failed: {}; reason: {}.\n\nContinue with next pipeline ..."
                msg = msg.format(self.name, e)
                self.exit_status = "FAILED"
          	raise PipelineExc(msg)
            inst.logger.removeHandler(fh)      
        logger.removeHandler(fh)
