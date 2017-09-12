import logging

logger = logging.getLogger(__name__)


def log_settings(settings):
    logger.info('\nSETTINGS')
    line = "-"*50
    logger.info(line)
    for key in settings:
        logger.info(" - {:20}{:20}".format(key, settings[key]))
    logger.info(line)


class PipelineExc(Exception):
    """
    only a specific pipeline failed, still continue with next pipelines
    """
    pass


class Fatal(Exception):
    """
    exception in the dsim run, all pipelines abort
    """
    pass

