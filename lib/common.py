import logging

logger = logging.getLogger(__name__)


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

