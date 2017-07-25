import sim_logger
import argparse
from lib import sim_reads
from lib import create_truth_vcf_and_fastas
from lib import get_dataset_ref_info_from_db
from lib import get_target_bed_from_db
from lib import parse_config
from lib import validate_and_parse_args


logger = sim_logger.logging.getLogger(__name__)


############################################################
# main - either generate truth VCF or simulate reads

def main(settings):

    # populate settings
    get_dataset_ref_info_from_db(settings)
    get_target_bed_from_db(settings)

    logger.info('\nSETTINGS')
    line = "-"*50
    logger.info(line)
    for key in settings:
        logger.info(" - {:20}{:20}".format(key, settings[key]))
    logger.info(line)

    create_truth_vcf_and_fastas(settings)
    sim_reads(settings)


############################################################
if __name__ == '__main__':
    '''
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', action='store', required=True)
    parser.add_argument('--varrate', action='store',
                        default=None, required=False)

    # parse args
    args = parser.parse_args()
    settings = {}
    parse_config(settings)
    validate_and_parse_args(args, settings)
    main(settings)
