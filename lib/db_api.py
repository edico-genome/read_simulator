"""
Methods to GET/ POST from DB/ LOOKUP TABLE
"""
import requests
import logging
import json
import os


logger = logging.getLogger(__name__)


class DBAPI(object):

    # lookup table for shared results
    # must be available across modules
    outputs_dict = {}

    def __init__(self, dataset_name):
        self.dataset_name = dataset_name
        self.outputs_dict[dataset_name] = {}

    ###################################
    # GET METHODS

    @staticmethod
    def get_from_db(__dataset_name, key_name):
        url = "http://data.edicogenome.com/api/get"
        __filter = {'name': __dataset_name, 'get_x': key_name}
        res = requests.get(url, params=__filter)
        res.raise_for_status()
        return res.text

    def get_dataset_ref_info(self):
        try:
            ref_type = self.get_from_db(self.dataset_name, 'ref_type')
            # use ref type as ref dataset name
            fasta = self.get_from_db(ref_type, 'fasta_file')
            dbsnp = self.get_from_db(ref_type, 'dbsnp_file')
            return {"ref_type": ref_type, "fasta": fasta, "dbsnp": dbsnp}
        except:
            logger.warning('failed to extract ref info from db')
            raise

    def get_target_bed(self, required=False):
        bed = self.get_from_db(self.dataset_name, 'target_region_bed')
        if bed:
            logger.info("Detected target regions bed: {}".format(bed))
        else:
            logger.warning("No target regions bed detected: {}".format(bed))
            if required:
                sys.exit(1)
        return bed

    def get_fastas(self):
        fasta0 = self.outputs_dict[self.dataset_name].get('fasta0')
        fasta1 = self.outputs_dict[self.dataset_name].get('fasta1')
        for f in [fasta0, fasta1]:
            if f:
                logger.info("Detected fasta: {}".format(f))
            else:
                logger.info("Did not detect fasta: {}".format(f))
        return {'fasta0': fasta0, 'fasta1': fasta1}

    ###################################
    # UPLOAD METHODS
    def set_fastas(self, fasta0, fasta1):
        for f in [fasta0, fasta1]:
            if not os.path.isfile(fasta0):
                raise Exception("Not a valid fasta: {}".format(f))
        self.outputs_dict[self.dataset_name] = {}
        self.outputs_dict[self.dataset_name]['fasta0'] = fasta0
        self.outputs_dict[self.dataset_name]['fasta1'] = fasta1

    def upload_to_db(self, key_name, value):
        url = "http://data.edicogenome.com/api/set"
        data = {
            'set_x': key_name,
            'name': self.dataset_name,
            'value': value
            }

        logger.info("Uploading to db. \nUrl: {}, \nData: {}".format(url, data))
        res = requests.post(url, params=data)
        res.raise_for_status()
        logger.info("Upload complete")

    def post_reads(self, fq1_path, fq2_path):
        self.upload_to_db('fastq_location_1', fq1_path)
        self.upload_to_db('fastq_location_2', fq2_path)

    def post_truth_vcf(self, truth_vcf):
        self.upload_to_db('truth_set_vcf', truth_vcf)
        self.upload_to_db('gatk_vcf_pre_filtered', truth_vcf)


    @staticmethod
    def post_requests(data):
        """
        may be related to SAVED OUTPUTS, or to CREATE DATASET
        :return:
        """
        headers = {'content-type': 'application/json'}
        r = requests.post(
            "http://data.edicogenome.com/api/dataset/submit",
            data=json.dumps(data),
            headers=headers
            )
        logging.info("status code: {}".format(r.status_code))
        logging.info("reason: {}".format(r.reason))
        r.raise_for_status()
