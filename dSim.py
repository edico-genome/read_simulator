#!/usr/bin/env python

import venv_check
import sys
import argparse
from lib import sim_logger
from lib.settings import Settings
from pipelines import pipelines


logger = sim_logger.logging.getLogger(__name__)


def is_venv():
    """ confirm we are running in a virtualenv """
    return (hasattr(sys, 'real_prefix') or
            (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix))


def pipeline_factory(pipeline_settings):
    """
    use pipeline settings ( pipeline name ) to determine which pipeline to instantiate
    """
    ThisPipelineClass = None
    ThisPipelineClass = getattr(pipelines, pipeline_settings["pipeline_name"])
    assert issubclass(ThisPipelineClass, pipelines.PipelinesBase), "This is not a valid pipeline: {}".format(c.name)


    if not ThisPipelineClass:
        logger.error("Please ensure pipeline: {} is registered and valid".format(
            pipeline_settings['pipeline_name']))
        sys.exit(1)

    p = ThisPipelineClass(pipeline_settings)
    return p


def instantiate_pipelines():
    """
    instantiate pipelines and validate pipeline settings
    """
    pipelines = []
    logger.info("\n\nValidating pipelines")
    for pipeline_settings in settings.runs:
        # default on
        if pipeline_settings.get("on/off", "on") in [0, None, "off"]:
            continue
        else:
            p = pipeline_factory(pipeline_settings)
            pipelines.append(p)
    return pipelines


def run_pipelines(pipelines):
    logger.info("\n\nRunning pipelines")
    for pipeline in pipelines:
        pipeline.run()


############################################################
# main
def main():
    pipelines = instantiate_pipelines()
    run_pipelines(pipelines)


############################################################
if __name__ == '__main__':
    """
    """
    if not is_venv():
        print "Please run in appropriate Python2.7 venv"
        sys.exit(1)

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run_settings', action='store', required=True)
    args = parser.parse_args()

    settings = Settings(args.run_settings)
    main()
    logger.info("Simulator Complete")
