#!./venv/bin/python

from lib import ensure_running_in_venv_upon_import
import sys, copy, argparse, os
from lib import sim_logger
from lib.settings import Settings
from pipelines import pipelines
from lib.common import PipelineExc, MPdb
from multiprocessing import Lock, Process, Queue
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


def instantiate_pipelines(settings):
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

def run_this_pipeline(p, result_queue=None):
    try:     
        p.run()
        state = "COMPLETED"
    except PipelineExc:
        msg = "Handled Exception: {}\nContinue with next pipeline ..."
        msg = msg.format(p.name)
        p.logger.error(msg, exc_info=True)        
        state = "HANLDED FAIL"
    except Exception:
        msg = "Unexpected Exception: {}\nContinue with next pipeline ..."
        msg = msg.format(p.name)
        p.logger.error(msg, exc_info=True)
        state = "UNEXPECTED FAIL"

    if result_queue:
        result_queue.put((p.name, p.dataset_name, state))

def run_pipelines(pipelines, run_serial):
    # run processes in parallel
    result_queue = Queue()
    MAX_PROCESSES = 2
    logger.info("\nRUNNING PIPELINES\n")    
    processes = []

    # run serial
    if run_serial:
        logger.info("Running in serial")
        for pipeline in pipelines:
            run_this_pipeline(pipeline)
        return 0

    # run parallel
    logger.info("Running in parallel")
    for pipeline in pipelines:
        p = Process(target=run_this_pipeline, args=(pipeline, result_queue))
        p.start()
        processes.append(p)
        if len(processes) == MAX_PROCESSES:
            for p in processes: p.join()
            processes = []
    if processes:
        for p in processes: p.join()

    logger.info("\nSIMULATOR PIPELINE EXIT STATE\n")
    exit_sig = 0
    for pipeline in pipelines:
        rv = result_queue.get()
        print("- {:10} {:35} {:15}".format(rv[0], rv[1], rv[2]))
        if rv[2] != "COMPLETED":
            exit_sig = 1
    return exit_sig
    

############################################################
# main
def main(settings, run_serial):
    pipelines_to_run = instantiate_pipelines(settings)
    return run_pipelines(pipelines_to_run, run_serial)
    

############################################################
if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run_settings', action='store', required=True)
    parser.add_argument('-s', '--run_serial', action='store_true', help='Default run parallel')
    args = parser.parse_args()

    logger.info("\nINPUT SETTINGS\n")
    settings = Settings(args.run_settings)
    settings.print_settings()
    main(settings, args.run_serial)

