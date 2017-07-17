import sim_logger
import pyaml
import argparse
from lib import create_truth_vcf
from lib import sim_reads
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

    logger.info('SETTINGS')
    line = "-"*50;
    logger.info("\n{}\n{}\n{}\n".format(line, pyaml.dump(settings), line))

    if settings['mode'] == 'truth_vcf':
        create_truth_vcf(settings)

    elif settings['mode'] == 'sim_reads':
        sim_reads(settings)

    else:
        create_truth_vcf(settings)
        sim_reads(settings)


############################################################
if __name__ == '__main__':
    '''
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', action='store', required=True)
    parser.add_argument('--varrate', action='store',
                        default=None, required=False)
    parser.add_argument('--mode', action='store', required=False,
                        default=None, help='Valid modes: truth_vcf/ sim_reads')

    # parse args
    args = parser.parse_args()
    settings = {}
    parse_config(settings)
    validate_and_parse_args(args, settings)
    main(settings)
