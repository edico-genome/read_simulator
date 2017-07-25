import os
import sys
import logging
import requests
import subprocess
import ConfigParser
import glob
from create_truth_vcf_and_modified_fastas import create_truth_vcf_and_modified_fastas

logger = logging.getLogger(__name__)


############################################################
# run setup functions

def parse_config(settings):
    config = ConfigParser.ConfigParser()
    config.read("config")
    settings['outdir'] = config.get('Paths', 'outdir')
    settings['PE100'] = config.get('Profiles', 'PE100')
    settings['indels'] = config.get('Profiles', 'indels')
    settings['gcdep'] = config.get('Profiles', 'gcdep')
    settings['varrate'] = config.get('Options', 'varrate')
    settings['errorrate'] = config.get('Options', 'errorrate')
    settings['useIndelErrors'] = config.get('Options', 'useIndelErrors')
    settings['reference'] = config.get('Options', 'reference')
    settings['reference'] = config.get('Options', 'reference')


def validate_and_parse_args(args, settings):
    settings['dataset'] = args.dataset

    if args.varrate:
        settings['varrate'] = args.varrate

    try:
        assert(float(settings['varrate']) <= 0.01)
    except:
        logger.error('Variant rate must be <= 0.01')
        sys.exit(1)

    run_name = "{}_indelErrors{}_errorRate{}_varRate{}".format(
        settings['dataset'], settings['useIndelErrors'],
        settings['errorrate'], settings['varrate'])
    settings['outdir'] = os.path.join(settings['outdir'], run_name)

    if not os.path.isdir(settings['outdir']):
        os.makedirs(settings['outdir'])


def get_dataset_ref_info_from_db(settings):
    try:
        settings['ref_type'] = get_from_db(settings['dataset'], 'ref_type')
        settings['fasta'] = get_from_db(settings['ref_type'], 'fasta_file')
        settings['dbsnp'] = get_from_db(settings['ref_type'], 'dbsnp_file')
    except:
        print 'failed to extract ref info from db'
        raise


def get_target_bed_from_db(settings):
    bed = get_from_db(settings['dataset'], 'target_region_bed')
    if bed:
        logger.info("Detected target regions bed: {}".format(bed))
        settings['target_bed'] = bed
    else:
        logger.info("No target regions bed detected: {}".format(bed))


############################################################
# common functions


def post_requests(data):
    headers = {'content-type': 'application/json'}
    r = requests.post(
        "http://data.edicogenome.com/api/dataset/submit",
        data=json.dumps(data),
        headers=headers)
    logging.info("status code: {}".format(r.status_code))
    logging.info("reason: {}".format(r.reason))
    r.raise_for_status()


def get_from_db(dataset_name, key_name):
    url = "http://data.edicogenome.com/api/get"
    filter = {'name': dataset_name,
              'get_x': key_name}
    res = requests.get(url, params=filter)
    res.raise_for_status()
    return res.text


def upload_to_db(key_name, dataset_name, value):
    url = "http://data.edicogenome.com/api/set"
    data = {
        'set_x': key_name,
        'name': dataset_name,
        'value': value
    }

    logging.info("Uploading to db. \nUrl: {}, \nData: {}".format(url, data))
    res = requests.post(url, params=data)
    try:
        res.raise_for_status()
        logging.info("Upload complete")
    except:
        logging.info("Upload failed")
        raise

############################################################
# vcf create


def bgzip_and_index_vcf(file_path):
    if '.gz' in file_path:
        raise Exception('Input VCF must be uncompressed')

    cmd = "bgzip -f " + file_path
    logger.info('bgzip cmd: {}'.format(cmd))
    subprocess.check_output(cmd, shell=True)

    vcf_gz = "{}.gz".format(file_path)
    cmd = 'tabix -f {}'.format(vcf_gz)
    logger.info('tabix cmd: {}'.format(cmd))
    subprocess.check_output(cmd, shell=True)
    return vcf_gz


def create_truth_vcf_and_fastas(settings):
    variant_settings = {}
    variant_settings['bed'] = settings['target_bed']
    variant_settings['outdir'] = settings['outdir']
    variant_settings['variant_rate'] = settings['varrate']
    variant_settings['ref_fasta'] = settings['fasta']


    truth_vcf, settings['modified_fasta_1'], settings['modified_fasta_2'] = \
        create_truth_vcf_and_modified_fastas(variant_settings)

    # compressed_vcf = bgzip_and_index_vcf(truth_vcf) #  this was only required for bcftools

    upload_to_db('truth_set_vcf', settings['dataset'], truth_vcf)
    upload_to_db('gatk_vcf_pre_filtered', settings['dataset'], truth_vcf)

    settings['truth_vcf'] = truth_vcf


############################################################
# simulate reads

def sim_reads(settings):
    logger.info('\nSIMULATE READS ...')
    # settings['truth_vcf'] = get_from_db(settings['dataset'], 'truth_set_vcf')
    run_pirs(settings)
    update_db_with_reads(settings)


def update_db_with_reads(settings):
    dataset = settings['dataset']
    upload_to_db('fastq_location_1', dataset, settings['fq1'])
    upload_to_db('fastq_location_2', dataset, settings['fq2'])


def run_pirs(settings):

    log = os.path.join(settings['outdir'], "pirs.log")

    options = {'PE100': settings['PE100'],
               'indels': settings['indels'],
               'gcdep': settings['gcdep'],
               'errorrate': settings['errorrate'],
               'f1': settings['modified_fasta_1'],
               'f2': settings['modified_fasta_2'],
               'outdir_pirs': os.path.join(settings['outdir'], 'pirs')}

    cmd = "pirs simulate -l 100 -x 30 -o {outdir_pirs}" + \
          " --insert-len-mean=180 --insert-len-sd=18 --diploid " + \
          " --base-calling-profile={PE100}" + \
          " --indel-error-profile={indels}" + \
          " --gc-bias-profile={gcdep}" + \
          " --phred-offset=33 --no-gc-bias -c gzip " + \
          " -t 48 {f1} {f2} "

    cmd = cmd.format(**options)

    if not settings['useIndelErrors']:
        cmd += " --no-indel-errors"
    cmd += " &> {}".format(log)

    logger.info("pirs cmd: {}".format(cmd))
    try:
        subprocess.check_output(cmd, shell=True)        
    except Exception() as e:
        logging.error('Error message %s' % e)
        raise

    logger.info('find output fastqs')
    fq1_list = glob.glob(os.path.join(settings['outdir'],'*1.fq.gz'))
    fq2_list = glob.glob(os.path.join(settings['outdir'],'*2.fq.gz'))
    settings['fq1'] = fq1_list[0]
    settings['fq2'] = fq2_list[0]
    update_db_with_reads(settings)
