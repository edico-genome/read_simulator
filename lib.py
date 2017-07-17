import os
import sys
import logging
import requests
import subprocess
import ConfigParser
import glob

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
    settings['mode'] = args.mode
    settings['dataset'] = args.dataset

    if args.varrate:
        settings['varrate'] = args.varrate

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
    logger.info('tabix cmd:     {}'.format(cmd))
    subprocess.check_output(cmd, shell=True)
    return vcf_gz


def intersect(settings):
    if settings['target_bed']:
        logger.info('Bedtools loading targed bed {}'
                    .format(settings['target_bed']))

        logger.info('Bedtools loading dbsnp {}'.format(settings['dbsnp']))

        intersect_vcf = os.path.join(settings['outdir'],
                                     'intersected_dbsnp.vcf')

        options = {'a': settings['dbsnp'],
                   'b': settings['target_bed'],
                   'vcf': intersect_vcf}

        cmd = "bedtools intersect -header -a {a} -b {b} > {vcf}".format(**options)

        logger.info('Intersect dbsnp with target_bed, writing result to: {}'
                    .format(intersect_vcf))
        logger.info('Bedtools cmd: {}'.format(cmd))

        subprocess.check_output(cmd, shell=True)
        settings['dbsnp'] = intersect_vcf


def create_truth_vcf(settings):

    logger.info('GENERATE TRUTH VCF')
    intersect(settings)

    desired_interval = int(1/float(settings['varrate']))
    logger.info('Estimated interval between variants: {}'
                .format(desired_interval))

    downsampled_vcf = os.path.join(
        settings['outdir'],
        'intersected_varrate_' + settings['varrate'] + '.vcf')

    logger.info('Writing downsampled VCF to: {}'
                .format(downsampled_vcf))

    previous_chromosome = None
    this_variant_position = None

    if '.gz' in settings['dbsnp']:
        import gzip
        open_gz_safe = gzip.open
    else:
        open_gz_safe = open

    with open_gz_safe(settings['dbsnp']) as stream_in, open(
            downsampled_vcf, 'w') as stream_out:
        for i, line in enumerate(stream_in):
            # skip headers
            if line[0] == '#':
                stream_out.write(line)
                continue

            this_chromosome = line.split()[0]
            this_variant_position = int(line.split()[1])

            if this_chromosome != previous_chromosome:
                previous_chromosome = this_chromosome
                ref_position = this_variant_position
                stream_out.write(line)
                belowline = line
                below = this_variant_position
                target = ref_position + desired_interval
                continue

            if this_variant_position > target:
                above = this_variant_position
                if (above - target) < (target - below):
                    stream_out.write(line)
                    # print("interval a: ", above - ref_position)
                    # print(line)
                    ref_position = above
                else:
                    # do not print same line more than once
                    if below > ref_position:
                        stream_out.write(belowline)
                        # print("interval b: ", below - ref_position)
                        # print(line)
                        ref_position = below
                    # take top line regardless
                    else:
                        stream_out.write(line)
                        # print("interval a reg: ", above - ref_position)
                        # print(line)
                        ref_position = above
                target = ref_position + desired_interval
            else:
                below = this_variant_position
                belowline = line

    logging.info("Created truth VCF, need to compress and create index")
    compressed_vcf = bgzip_and_index_vcf(downsampled_vcf)
    upload_to_db('truth_set_vcf', settings['dataset'], compressed_vcf)
    settings['truth_vcf'] = compressed_vcf


############################################################
# simulate reads

def sim_reads(settings):
    logger.info('SIMULATE READS')
    settings['truth_vcf'] = get_from_db(settings['dataset'], 'truth_set_vcf')
    generate_fasta(settings, 1)
    generate_fasta(settings, 2)
    run_pirs(settings)
    update_db_with_reads(settings)


def update_db_with_reads(settings):
    dataset = settings['dataset']
    outdir = settings['outdir']
    upload_to_db('fastq_location_1', dataset,
                 settings['fq1'])
    upload_to_db('fastq_location_2', dataset,
                 settings['fq2'])


def generate_fasta(settings, haplotype):

    # index = out_dir_root + 'mutations.vcf.gz.gzi'
    newfasta = '{}/modified_{}.fa'.format(settings['outdir'], haplotype)
    liftover = os.path.join(settings['outdir'],
                            "liftover_{}.txt".format(haplotype))
    varerror = os.path.join(settings['outdir'],
                            "varerror_{}.txt".format(haplotype))

    options = {'liftover': liftover,
               'hap': haplotype,
               'fasta': settings['fasta'],
               'vcf': settings['truth_vcf'],
               'new_fasta': newfasta,
               'varerror': varerror}

    cmd = "bcftools consensus -c {liftover} " + \
          "-H {hap} -f {fasta} " + \
          "{vcf} 1> {new_fasta} 2> {varerror}"

    cmd = cmd.format(**options)

    logger.info("Generate fasta for haplotype {}".format(haplotype))
    logger.info("bcftools cmd: {}".format(cmd))

    try:
        subprocess.check_output(cmd, shell=True)
        settings['modified_fasta_{}'.format(haplotype)] = newfasta
    except Exception() as e:
        logging.error('Failed to run bcftools %s' % e)
        raise


def run_pirs(settings):

    log = os.path.join(settings['outdir'], "pirs.log")

    options = {'PE100': settings['PE100'],
               'indels': settings['indels'],
               'gcdep': settings['gcdep'],
               'errorrate': settings['errorrate'],
               'f1': settings['modified_fasta_1'],
               'f2': settings['modified_fasta_2'],
               'outdir_pirs': os.path.join(settings['outdir'], 'pirs')}

    # " --error-rate={errorrate}" 
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
