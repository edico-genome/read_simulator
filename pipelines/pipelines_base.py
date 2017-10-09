import os
from abc import ABCMeta, abstractproperty
from lib.db_api import DBAPI
from lib.common import PipelineExc
import logging, sys

logger = logging.getLogger(__name__)


###########################################################
# all pipelines must inherit from this base

class PipelinesBase(object):
    """
    Base class for all other classes in this module
    """
    __metaclass__ = ABCMeta

    def __init__(self, settings):
        msg = "validating pipeline: {}".format(self.__class__.__name__)
        logger.info(msg)
        self.exit_status = "FAILED"
        self.pipeline_settings = settings
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
        fh = logging.FileHandler(os.path.join(self.pipeline_settings['outdir'],
                                              self.pipeline_settings['pipeline_name']))
        logger.addHandler(fh)
        logger.info("\nPIPELINE: {}".format(self.name))
        for inst in self.module_instances:
            try:
                inst.before_run()
                inst.run()
                inst.after_run()
                self.exit_status = "COMPLETED"
            except PipelineExc as e:
                logger.error("Fatal error: {}".format(e), exc_info=True)
                msg = "Pipeline Failed: {}; \nContinue with next pipeline ..."
                msg = msg.format(self.name)
                PipelineExc(msg)
                return
            except Exception as e:
                logger.error("Fatal error: {}".format(e), exc_info=True)
                sys.exit(1)
        logger.removeHandler(fh)
