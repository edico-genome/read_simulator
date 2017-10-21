#!/usr/bin/env python

from lib import ensure_running_in_venv_upon_import
import sys, copy, argparse, os
from lib import sim_logger
from lib.settings import Settings
from pipelines import pipelines
from lib.common import PipelineExc
from multiprocessing import Lock, Process
import logging

# logger will be used for simulator level logging
# logging will be used for pipeline level logging
logger = sim_logger.logging.getLogger(__name__)


############################################################
# main
def pipeline_factory(pipeline_name):
    """
    use pipeline settings ( pipeline name ) to determine which pipeline to instantiate
    """
    ThisPipelineClass = getattr(pipelines, pipeline_name)

    if not ThisPipelineClass and issubclass(ThisPipelineClass, pipelines.PipelinesBase):
        logger.error("Invalid pipeline: {}".format(pipeline_name))
        sys.exit(1)

    return ThisPipelineClass


def instantiate_pipelines():
    """
    instantiate pipelines and validate pipeline settings
    """
    pipelines = []
    # lock to manage race parallel processes race conditions 
    lock = Lock()

    logger.info("\nVALIDATING PIPELINES\n")
    for p_idx, pipeline_settings in enumerate(settings.runs):

        # turn a pipeline off by specifying num_runs as 0
        num_runs = pipeline_settings.get("num_runs", 0)

        # start_idx determines the first dataset name's starting idx
        start_idx = pipeline_settings.get("start_idx", 0)

        if num_runs:
            logger.info("Validating run: {}\n".format(p_idx))
        else:
            logger.info("Skipping run: {}\n".format(p_idx))
            
        for idx in range(start_idx, start_idx + num_runs):           
            logger.info("Pipeline sub index: {}\n".format(idx))
            # class factory and instantiate pipeline object
            Pipeline = pipeline_factory(pipeline_settings["pipeline_name"])
            p = Pipeline(pipeline_settings, idx)
            
            # give each pipeline an idependent logger
            log_name = "dSim_{}".format(p.pipeline_settings["dataset_name"])
            log_path = os.path.join(p.pipeline_settings["outdir"],
                                    p.pipeline_settings["dataset_name"]+'.log')
            fh = logging.FileHandler(log_path, mode='w')
            fh.setLevel(logging.DEBUG)
            format = "%(asctime)-6s: %(name)s - %(levelname)s - %(message)s"
            fmt = logging.Formatter(format)
            fh.setFormatter(fmt)
            local_logger = logging.getLogger(log_name)
            local_logger.addHandler(fh)
            logger.info("Init local logging: {}".format(log_path))
            p.logger = local_logger

            # pipeline/ dataset directory
            p.pipeline_settings["lock"] = lock

            # validate all submodules for each pipeline is ready (use local logger) 
            p.instantiate_modules()

            # append to list of instantiated pipelines
            pipelines.append(p)
    return pipelines


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

    # run processes in parallel
    MAX_PROCESSES = 5
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


def print_summary(pipelines):
    # print summary
    logger.info("\nSIMULATOR SUMMARY\n")
    for pipeline in pipelines:
        logger.info("{} {} {}".format(
                pipeline.name,
                pipeline.dataset_name,
                pipeline.exit_status))
        if pipeline.exit_status != "COMPLETED":
            for m in pipeline.module_instances:
                logger.info(" - {} {}".format(
                        m.name,
                        pipeline.exit_status))
        

############################################################
# main
def main():
    pipelines_to_run = instantiate_pipelines()
    run_pipelines(pipelines_to_run)
    print_summary(pipelines_to_run)


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

