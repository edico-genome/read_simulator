#!/usr/bin/env python

import sys, copy
import argparse
from lib import sim_logger
from lib.settings import Settings
from pipelines import pipelines


logger = sim_logger.logging.getLogger(__name__)


def ensure_running_in_venv():
    """ confirm that we are running in a virtualenv """
    test1 = hasattr(sys, 'real_prefix')
    test2 = hasattr(sys, 'base_prefix') and (sys.base_prefix != sys.prefix)
    if not (test1 or test2):
        print "Please run in appropriate Python2.7 venv"
        sys.exit(1)


def pipeline_factory(pipeline_name):
    """
    use pipeline settings ( pipeline name ) to determine which pipeline to instantiate
    """
    ThisPipelineClass = getattr(pipelines, pipeline_name)

    if not ThisPipelineClass:
        logger.error("Please ensure pipeline: {} is registered"
                     .format(pipeline_name))
        sys.exit(1)

    assert issubclass(ThisPipelineClass, pipelines.PipelinesBase),\
        "Please ensure pipeline: {} is a valid subclass of type pipeline: {}"\
        .format(pipeline_name)

    return ThisPipelineClass


def instantiate_pipelines():
    """
    instantiate pipelines and validate pipeline settings
    """
    pipelines = []
    logger.info("\nVALIDATING PIPELINES\n")
    for pipeline_settings in settings.runs:
        if pipeline_settings.get("on/off", "on") in [0, None, "off"]:
            continue
        else:
            num_runs = pipeline_settings.get("num_runs", "1")
            for idx in range(num_runs):
                this_pipeline_settings = copy.deepcopy(pipeline_settings)
                if idx > 0:
                    this_pipeline_settings["dataset_name"] += "_{}".format(idx)
                Pipeline = pipeline_factory(this_pipeline_settings["pipeline_name"])
                p = Pipeline(this_pipeline_settings)
                pipelines.append(p)
    return pipelines


def run_pipelines(pipelines):
    logger.info("\nRUNNING PIPELINES\n")
    for pipeline in pipelines:
        pipeline.run()

    logger.info("\nSIMULATOR SUMMARY\n")
    for pipeline in pipelines:
        logger.info("{} {} {}".format(
                pipeline.name,
                pipeline.dataset_name,
                pipeline.exit_status))


############################################################
# main
def main():
    pipelines = instantiate_pipelines()
    run_pipelines(pipelines)


############################################################
if __name__ == '__main__':
    """
    """
    ensure_running_in_venv()
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run_settings', action='store', required=True)
    args = parser.parse_args()

    logger.info("\nINPUT SETTINGS\n")
    settings = Settings(args.run_settings)

    main()

