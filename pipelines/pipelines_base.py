from abc import ABCMeta, abstractproperty
from lib.db_api import DBAPI
import logging

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
        logger.info("-" * 60)
        logger.info(msg)
        self.pipeline_settings = settings
        self.module_instances = []
        self.instantiate_modules()
        self._class_name = self.__class__.__name__
        self.db_api = DBAPI(self.pipeline_settings["dataset_name"])
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
            inst = C(self.pipeline_settings)
            self.module_instances.append(inst)

    def run(self):
        msg = "running pipeline: {}".format(self.name)
        logger.info("-" * 60)
        logger.info(msg)
        for inst in self.module_instances:
            inst.run()

