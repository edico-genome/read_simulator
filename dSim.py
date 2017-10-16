#!/usr/bin/env python

import os
from lib import ensure_running_in_venv_upon_import
import sys, copy, argparse
from lib import sim_logger
from lib.settings import Settings
from pipelines import pipelines
from lib.common import PipelineExc
from multiprocessing import Lock, Process

logger = sim_logger.logging.getLogger(__name__)

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
    _pipelines = []
    lock = Lock()

    logger.info("\nVALIDATING PIPELINES\n")
    for pipeline_settings in settings.runs:
        if pipeline_settings.get("on/off", "off") in [0, None, "off"]:
            continue
        else:
            num_runs = pipeline_settings.get("num_runs", 1)
            for idx in range(num_runs):
                this_pipeline_settings = copy.deepcopy(pipeline_settings)
                if idx > 0:
                    this_pipeline_settings["dataset_name"] += "_{}".format(idx)
                this_pipeline_settings["workdir"] += "/p_{}".format(idx)
                if not os.path.isdir(this_pipeline_settings["workdir"]):
                    os.makedirs(this_pipeline_settings["workdir"])
                this_pipeline_settings["lock"] = lock
                Pipeline = pipeline_factory(this_pipeline_settings["pipeline_name"])
                p = Pipeline(this_pipeline_settings)
                _pipelines.append(p)
    return _pipelines


def run_this_pipeline(_pipeline):
    try:
        _pipeline.run()
    except PipelineExc as e:
        logger.error("Pipeline failed: {}".format(e), exc_info=True) 
        raise
    except Exception as e:
        logger.error("Unknown fatal error: {}".format(e), exc_info=True)
        raise


def run_pipelines(pipelines):
    MAX_PROCESSES = 10

    logger.info("\nRUNNING PIPELINES\n")    
    processes = []
    for pipeline in pipelines:
        p = Process(target=run_this_pipeline, args=(pipeline,))
        p.start()
        processes.append(p)
        if len(processes) == MAX_PROCESSES:
            for p in processes: p.join()
            processes = []
    for p in processes: p.join()

    logger.info("\nSIMULATOR SUMMARY\n")
    for pipeline in pipelines:
        logger.info("{} {} {}".format(
                pipeline.name,
                pipeline.dataset_name,
                pipeline.exit_status))


############################################################
# main
def main():
    pipelines_to_run = instantiate_pipelines()
    run_pipelines(pipelines_to_run)


############################################################
if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run_settings', action='store', required=True)
    args = parser.parse_args()

    logger.info("\nINPUT SETTINGS\n")
    settings = Settings(args.run_settings)
    settings.print_settings()
    main()

