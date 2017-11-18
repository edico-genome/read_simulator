import os
import sys
import gzip
import logging
import random
import subprocess
from copy import deepcopy
from collections import OrderedDict


logger = logging.getLogger(__name__)


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


def write_vcf_header(stream):
    stream.write('##fileformat=VCFv4.2\n')
    stream.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')


def write_vcf_line(variant, stream_out):
    vcf_line = "{chr}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t." + \
               "\t{info}\t{format}\t{genotype}\n"
    vcf_line = vcf_line.format(**variant)
    stream_out.write(vcf_line)


def print_vcf(settings):
    """ create normalized truth vcf """

    truth_vcf = os.path.join(settings['outdir'], 'truth.vcf')
    norm_split_vcf = os.path.join(settings['outdir'], 'norm_split_truth.vcf')
    norm_merged_vcf_gz = os.path.join(settings['outdir'], 'norm_merged_truth.vcf.gz')

    logger.info('Generate truth VCF: {}'.format(truth_vcf))
    with open(truth_vcf, 'w') as stream_out:
        write_vcf_header(stream_out)

        for chr in settings['parsed_fasta']:
            print("Adding chr {} to VCF".format(chr))
            if chr not in settings['var_info_for_vcf']:
                print("- no variants in chr {}- continue".format(chr))
                continue

            # loop over all variants in chr
            for var_info in reversed(settings['var_info_for_vcf'][chr]):
                genotype = [None, None]
                pos, ref_0, ref_1, alt_0, alt_1 = var_info
                assert ((ref_0[0] == ref_1[0]))

                ref_is_lower_case= False
                if ref_0[0].lower() not in ['a', 'b', 'c', 'd']:
                    continue

                if ref_0[0] in ['a', 'b', 'c', 'd']:
                    ref_is_lower_case = True

                if ref_is_lower_case:
                    alt_0 = alt_0.lower()
                    alt_1 = alt_1.lower()
                    
                # variant template
                qual = '.'
                info = '.'
                format = 'GT:AD:DP:GQ:PL:SB'
                genotype_template = '{}:1,1000:1000:10:10000,1000,0:0,1,1000,1000'

                # only print actual variants
                if ((alt_0 == alt_1) and (ref_0 == alt_0) and (ref_1 == ref_0)):
                    continue

                # normalize variants ( e.g. for heter insertion and deletion ) 
                geno = "unknown" 
                alleles = []
                longest_ref = ref_0 if len(ref_0) > len(ref_1) else ref_1

                for ref, alt in zip((ref_0, ref_1), (alt_0, alt_1)):
                    if ref != longest_ref:
                        alt = alt + longest_ref[len(ref):]

                    if alt == longest_ref: # if this is an actual variant then the other allele must differ
                        geno = "0/1"
                    else:
                        if alt in alleles: # this allele was already defined - i.e. homozygous
                            geno = "1/1"
                        else:
                            alleles.append(alt)  
                            if len(alleles) == 2: # is this the second allele we're appending?
                                geno = "1/2"

                variant = {
                    'chr': chr,  
                    'pos': pos + 1,
                    'ref': longest_ref,
                    'alt': ",".join(alleles),
                    'qual': qual,
                    'format': format,
                    'info': info,
                    'genotype': genotype_template.format(geno)
                }
                write_vcf_line(variant, stream_out)
                
    # normalize split 
    cmd = "bcftools norm -f {} {} > {}".format(settings['fasta_file'], truth_vcf, norm_split_vcf)
    logger.info('{}'.format(cmd))
    subprocess.check_output(cmd, shell=True)

    # normalize merge
    cmd = "bcftools norm -m +any -N {} -O z > {}".format(norm_split_vcf, norm_merged_vcf_gz)
    logger.info('{}'.format(cmd))
    subprocess.check_output(cmd, shell=True)

    # set truth
    settings['truth_vcf'] = norm_merged_vcf_gz

    # index
    cmd = "bcftools index -f {}".format(norm_merged_vcf_gz)
    subprocess.check_output(cmd, shell=True)


def open_gz_safe(file_path):
    if '.gz' in file_path:
        return gzip.open(file_path)
    else:
        return open(file_path)


def sample_variant():
    vars = (
        ('snp', 'G'),
        ('snp', 'C'),
        ('snp', 'A'),
        ('snp', 'T')
        )
    r = random.randint(0, len(vars)-1)
    return vars[r]


def define_variants(settings):
    '''
    variants will be created in the specified bed regions
    variants will be created at specified intervals
    '''
    print('Sampling variants')

    # regions = ['chr1 0 50', 'chr1 100 150']
    bases_between_variants = int(1 / float(settings['varrate']))
    print('bases between variants: {}'.format(bases_between_variants))

    settings['sampled_vars'] = {}
    with open(settings['sorted_bed']) as regions:
        for reg in regions:
            _chr, _from, _to = reg.split()[0:3]
            _from = int(_from) + random.randint(0, 15)
            _to = int(_to)
            # avoid boundary effects
            if _to - 100 <= _from:
                logger.error('Please make sure regions are > 100bp')
                sys.exit(1)
            _from += 35 
            _to -= 25
            _pos = range(_from, _to, bases_between_variants)
            _allele1 = [sample_variant() for i in range(len(_pos))]
            _allele2 = [sample_variant() for i in range(len(_pos))]
            # _allele2 = [None for i in range(len(_pos))] # for heter allelles

            if _chr not in settings['sampled_vars']:
                settings['sampled_vars'][_chr] = []

            for _p, _a1, _a2 in zip(_pos, _allele1, _allele2):
                settings['sampled_vars'][_chr].append((_p, _a1, _a2))

    for chr in settings['sampled_vars']:
        logger.info('chr: {}, nr of variants {}'.format(chr, len(settings['sampled_vars'][chr])))


##################################################
def get_reference_allele_from_fasta(settings):

    for _chr in settings['sampled_vars']:
        # read fasta in reverse
        for var_info in reversed(settings['sampled_vars'][_chr]):

            index_in_chr_where_to_inject = var_info[0]
            alleles_0_1 = var_info[1:]

            # how to access the variant from global index
            fasta_line_index = int(
                index_in_chr_where_to_inject /
                settings['length_of_fasta_line'])
            fasta_line_rel_pos = index_in_chr_where_to_inject % \
                settings['length_of_fasta_line']

            alt = [None, None]
            ref = [None, None]

            this_fasta_line = settings['parsed_fasta'][_chr][fasta_line_index]
            try:
                next_fasta_line = settings['parsed_fasta'][_chr][fasta_line_index+1]
            except:
                next_fasta_line = None

            for haplotype in [0, 1]:
                # allele index 0 = type, index 1 = attribute
                allele = alleles_0_1[haplotype]

                # set this allele to ref
                if not allele:
                    ref[haplotype] = this_fasta_line[fasta_line_rel_pos]
                    alt[haplotype] = this_fasta_line[fasta_line_rel_pos]
                    continue

                # SNPs
                if allele[0] == 'snp':
                    #  get info for vcf
                    ref[haplotype] = this_fasta_line[fasta_line_rel_pos]
                    alt[haplotype] = allele[1]

                # INSERTs
                if allele[0] == 'ins':
                    #  get info for vcf
                    ref[haplotype] = this_fasta_line[fasta_line_rel_pos]
                    alt[haplotype] = ref[haplotype] + allele[1]

                # DELs
                if allele[0] == 'del':
                    # check if this deletion flows into the next line
                    deletion_end_index = fasta_line_rel_pos + allele[1]
                    len_this_line = len(this_fasta_line)
                    del_overflow_into_next_line = deletion_end_index - len_this_line

                    # if this deletion does not flow into next line
                    if del_overflow_into_next_line <= 0:
                        #  get info for vcf
                        # allele[1] indicates how long the deletion is
                        ref[haplotype] = this_fasta_line[fasta_line_rel_pos:fasta_line_rel_pos+allele[1]+1]
                        alt[haplotype] = ref[haplotype][0]
                    # if this deletion overflows into next line
                    else:
                        #  get info for vcf
                        # allele[1] indicates how long the deletion is
                        ref[haplotype] = this_fasta_line[fasta_line_rel_pos:fasta_line_rel_pos+allele[1]+1] + \
                                         next_fasta_line[:del_overflow_into_next_line+1]
                        alt[haplotype] = ref[haplotype][0]

            # store the information we'll need for the VCF
            if _chr not in settings['var_info_for_vcf']:
                settings['var_info_for_vcf'][_chr] = []
            settings['var_info_for_vcf'][_chr].append((
                index_in_chr_where_to_inject, ref[0], ref[1], alt[0], alt[1]))


##################################################
def parse_ref_fasta(settings):
    logger.info('parse reference fasta')

    with open_gz_safe(settings['fasta_file']) as stream:
        settings['length_of_fasta_line'] = None
        this_header = None
        settings['parsed_fasta'] = OrderedDict()
        this_chr_fasta = []
        counter = 0

        for line in stream:
            counter += 1
            line = line.replace('\n', '')

            if ((len(line.split()[0]) < 40) and ('>' in line)):
                # dump old chr
                if this_header:
                    settings['parsed_fasta'][this_header] = this_chr_fasta
                # transition to new chr
                this_header = line.split()[0].strip().replace(">", "")
                print('Parsing fasta chr: {}'.format(this_header))
                this_chr_fasta = []
                
            else:
                this_chr_fasta.append(line)
                if not settings['length_of_fasta_line']:
                    settings['length_of_fasta_line'] = len(line)
        # capture output from last chr
        settings['parsed_fasta'][this_header] = this_chr_fasta


##################################################
def create_truth_vcf(settings):
    """
    takes as input a FASTA file and a specification regarding which
    regions/ type and frequency of variants you want to simulate
    """

    settings['var_info_for_vcf'] = {}
    settings['mod_fasta'] = []

    logger.info('\nVLRD SIMULATING VARIANTS ...')
    parse_ref_fasta(settings)
    define_variants(settings)
    get_reference_allele_from_fasta(settings)
    print_vcf(settings)
    logger.info('variant simulation complete')

    # delete fasta from memory
    settings.pop('sampled_vars')

    # convert to consistent module names
    settings['truth_set_vcf'] = settings['truth_vcf']



"""
('snp','G'),
('snp','C'),
('snp','A'),
('snp','T'),
('snp','G'),
('snp','C'),
('snp','A'),
('snp','T'),
('snp','G'),
('snp','C'),
('snp','A'),
('snp','T'),
('snp','G'),
('snp','C'),
('snp','A'),
('snp','T'),
('snp','G'),
('snp','C'),
('snp','A'),
('snp','T'),
('snp','G'),
('snp','C'),
('snp','A'),
('snp','T'),
('snp','G'),
('snp','C'),
('snp','A'),
('snp','T'),
('ins','AT'),
('ins','GA'),
('ins','TCAT'),
('ins','TTAT'),
('ins','TAT'),
('ins','TTCAT'),
('ins','GTCCAT'),
('ins','AATCGCTGACGCAT'),
#('ins','CATTGGTCGGGGATTTCTGA'),
#('ins','TCTGCCGATTCATTGGTGCAGACTATCGGGGATTTCTGATGCA'),
('del',1),
('del',2),
('del',3),
('del',4),
('del',5),
('del',6),
('del',10),
#('del',20),
#('del',50),
"""
