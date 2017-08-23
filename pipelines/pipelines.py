"""
Classes for stitching modules into a pipeline
"""
from abc import ABCMeta, abstractproperty
from modules.modules import VLRDmod, Pirs, CompositeDataset
import logging

logger = logging.getLogger(__name__)


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

    @abstractproperty
    def modules(self):
        """
        :rtype: list
        """
        pass

    def instantiate_modules(self):
        """
        responsible for validating settings
        """
        for C in self.modules:
            logger.info("- validating module: {}".format(C.name))
            inst = C(self.pipeline_settings)
            self.module_instances.append(inst)

    def run(self):
        msg = "running pipeline: {}".format(self.__class__.__name__)
        logger.info("-" * 60)
        logger.info(msg)
        for inst in self.module_instances:
            inst.run()


###########################################################
# pipelines to create a single sample

class VLRDPipeline(PipelinesBase):
    modules = [VLRDmod, Pirs]


###########################################################
# pipelines to create composite samples

class CompositeDatasetPipeline(PipelinesBase):
    modules = [CompositeDataset]


###########################################################
# a dictionary of listed classes
"""
The run config should refer to pipelines using the keys in registered_pipelines
"""

registered_pipelines = {
    "vlrd": VLRDPipeline,
    "composite_dataset": CompositeDatasetPipeline
    }
