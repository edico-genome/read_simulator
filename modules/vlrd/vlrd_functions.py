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
    truth_vcf = os.path.join(settings['outdir'], 'truth.vcf')
    settings['truth_vcf'] = truth_vcf
    logger.info('Generate truth VCF: {}'.format(truth_vcf))

    with open(truth_vcf, 'w') as stream_out:
        write_vcf_header(stream_out)

        chrs = [c.replace('>','') for c in settings['parsed_fasta']]
        for chr in chrs:
            print("Adding chr {} to VCF".format(chr))
            if chr not in settings['var_info_for_vcf']:
                print("- continue")
                continue
            for var_info in reversed(settings['var_info_for_vcf'][chr]):
                genotype = [None, None]
                pos, ref_0, ref_1, alt_0, alt_1 = var_info
                assert ref_0[0] == ref_1[0]

                # variant template
                qual = '.'
                info = 'N/A'
                format = 'GT:AD:DP:GQ:PL:SB'
                genotype_template = '{}:1,1000:1000:10:10000,1000,0:0,1,1000,1000'

                # only print actual variants
                if ((alt_0 == alt_1) and (ref_0 == alt_0)):
                    continue

                # snps and inserts are easy
                alt = []
                genotype = [0, 0]
                if ((len(ref_0) == 1) and (len(ref_1) == 1)):
                    for idx, i in enumerate([alt_0, alt_1]):
                        if ref_0 != i:
                            alt.append(i)
                            genotype[idx] = len(alt)

                    variant = {
                        'chr': chr,  # offset 1
                        'pos': pos + 1,
                        'ref': ref_0,
                        'alt': ','.join(alt),
                        'qual': qual,
                        'format': format,
                        'info': info,
                        'genotype': genotype_template.format('|'.join([str(g) for g in genotype]))
                    }
                    write_vcf_line(variant, stream_out)

                # deletions - split over 2 lines
                if not ((len(ref_0) == 1) and (len(ref_1) == 1)):
                    for ref, alt, geno in zip((ref_0, ref_1), (alt_0, alt_1), ('0|1', '1|0')):
                        variant = {
                            'chr': chr,  # offset 1
                            'pos': pos + 1,
                            'ref': ref,
                            'alt': alt,
                            'qual': qual,
                            'format': format,
                            'info': info,
                        'genotype': genotype_template.format(geno)
                        }
                        write_vcf_line(variant, stream_out)


def open_gz_safe(file_path):
    if '.gz' in file_path:
        return gzip.open(file_path)
    else:
        return open(file_path)


def sample_variant():
    vars = (('snp', 'G'),
            ('snp', 'C'),
            ('snp', 'A'),
            ('snp', 'T'),
            ('ins', 'AT'),
            ('del', 2))

    # r = random.randint(0, len(vars)-1)
    r = random.randint(0, 4)
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
    with open(settings['target_bed']) as regions:
        for reg in regions:
            _chr, _from, _to = reg.split()[0:3]
            _from = int(_from) + random.randint(0, 20)
            _to = int(_to)
            # avoid boundary effects
            if _to - 100 <= _from:
                logger.error('Please make sure regions are > 100bp')
                sys.exit(1)
            _from += 50 
            _to -= 50
            _pos = range(_from, _to, bases_between_variants)
            _allele1 = [sample_variant() for i in range(len(_pos))]
            _allele2 = [sample_variant() for i in range(len(_pos))]

            if _chr not in settings['sampled_vars']:
                settings['sampled_vars'][_chr] = []

            for _p, _a1, _a2 in zip(_pos, _allele1, _allele2):
                settings['sampled_vars'][_chr].append((_p, _a1, _a2))

    for chr in settings['sampled_vars']:
        logger.info('chr: {}, nr of variants {}'.format(chr, len(settings['sampled_vars'][chr])))


def add_variants_to_fasta(settings):
    logger.info('adding variants to fasta')

    settings['mod_fasta'] = [None, None]
    print("Deep copy fasta 0")
    settings['mod_fasta'][0] = deepcopy(settings['parsed_fasta'])
    print("Deep copy fasta 1")
    settings['mod_fasta'][1] = deepcopy(settings['parsed_fasta'])

    for chr in settings['sampled_vars']:
        # update fasta in reverse so that indels
        # do not mix up the index positions for unmodified variants

        # fasta requires chr in '>'chr format
        f_chr = '>'+chr

        for var_info in reversed(settings['sampled_vars'][chr]):

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

            this_fasta_line = settings['parsed_fasta'][f_chr][fasta_line_index]
            try:
                next_fasta_line = settings['parsed_fasta'][f_chr][fasta_line_index+1]
            except:
                next_fasta_line = None

            for haplotype in [0, 1]:
                # allele index 0 = type, index 1 = attribute
                allele = alleles_0_1[haplotype]

                # SNPs
                if allele[0] == 'snp':
                    #  get info for vcf
                    ref[haplotype] = this_fasta_line[fasta_line_rel_pos]
                    alt[haplotype] = allele[1]
                    #  update fasta
                    settings['mod_fasta'][haplotype][f_chr][fasta_line_index] = \
                        this_fasta_line[:fasta_line_rel_pos] + \
                        allele[1] + \
                        this_fasta_line[fasta_line_rel_pos+1:]

                # INSERTs
                if allele[0] == 'ins':
                    #  get info for vcf
                    ref[haplotype] = this_fasta_line[fasta_line_rel_pos]
                    alt[haplotype] = ref[haplotype] + allele[1]
                    #  update fasta
                    settings['mod_fasta'][haplotype][f_chr][fasta_line_index] = \
                        this_fasta_line[:fasta_line_rel_pos+1] + \
                        allele[1] + \
                        this_fasta_line[fasta_line_rel_pos+1:]

                # DELs
                if allele[0] == 'del':

                    raise Exception("deletions not currently supported")

                    """
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

                        # update fasta
                        settings['mod_fasta'][haplotype][f_chr][fasta_line_index] = \
                            this_fasta_line[:fasta_line_rel_pos+1] + \
                            'D'*allele[1] + \
                            this_fasta_line[fasta_line_rel_pos+allele[1]+1:]

                    else:
                        #  get info for vcf
                        # allele[1] indicates how long the deletion is
                        ref[haplotype] = this_fasta_line[fasta_line_rel_pos:fasta_line_rel_pos+allele[1]+1] + \
                                         next_fasta_line[:del_overflow_into_next_line+1]
                        alt[haplotype] = ref[haplotype][0]

                        # first cut the deletion from the end of this line
                        settings['mod_fasta'][haplotype][f_chr][fasta_line_index] = \
                            this_fasta_line[:fasta_line_rel_pos+1] + \
                            'D'*(len_this_line - 1 - fasta_rel_index)

                        # then cut the overflow from the start of the next line
                        settings['mod_fasta'][haplotype][f_chr][fasta_line_index+1][fasta_line_rel_pos+1:] = \
                            'D'*del_overflow_into_next_line + \
                            next_fasta_line[del_overflow_into_next_line:]
                    """

            # provide the information we'll need for the VCF
            if chr not in settings['var_info_for_vcf']:
                settings['var_info_for_vcf'][chr] = []
            settings['var_info_for_vcf'][chr].append((
                index_in_chr_where_to_inject, ref[0], ref[1], alt[0], alt[1]))


def write_fasta_to_file(settings, haplotype):

    fasta_out = os.path.join(
        settings['outdir'], 'haplotype_{}.fasta'.format(haplotype))

    settings['mod_fasta_path_{}'.format(haplotype)] = fasta_out

    logger.info('write fasta {} to file {}'.format(haplotype, fasta_out))
    with open(fasta_out, 'w') as stream_out:
        # for chr in settings['mod_fasta'][haplotype]:
        for chr in settings['parsed_fasta']: #  ordered dict - preserve fasta order
            print("Writing chr {}".format(chr))
            # make everything of same line length
            stream_out.write("{}\n".format(chr))
            # import pdb; pdb.set_trace()
            for line in settings['mod_fasta'][haplotype][chr]:
                stream_out.write("{}\n".format(line))


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
            # if counter > 100000:
            #    print 'BREAK'
            #    break

            if ((len(line.split()[0]) < 40) and ('>' in line)):
                # dump old chr
                if this_header:
                    settings['parsed_fasta'][this_header] = this_chr_fasta
                # transition to new chr
                print('Parsing fasta chr: {}'.format(line))
                this_header = line.split()[0].strip()
                this_chr_fasta = []
            else:
                this_chr_fasta.append(line)
                if not settings['length_of_fasta_line']:
                    settings['length_of_fasta_line'] = len(line)
        # capture output from last chr
        settings['parsed_fasta'][this_header] = this_chr_fasta


def create_truth_vcf_and_fastas(settings):
    """
    takes as input a FASTA file and a specification regarding which
    regions/ type and frequency of variants you want to simulate
    """

    settings['var_info_for_vcf'] = {}
    settings['mod_fasta'] = []

    logger.info('\nVLRD SIMULATING VARIANTS ...')
    parse_ref_fasta(settings)
    define_variants(settings)
    add_variants_to_fasta(settings)
    print_vcf(settings)
    write_fasta_to_file(settings, 0)
    write_fasta_to_file(settings, 1)
    logger.info('variant simulation complete')

    res = {"truth_vcf": settings['truth_vcf'],
           "fasta0": settings['mod_fasta_path_0'],
           "fasta1": settings['mod_fasta_path_1']}

    return res



