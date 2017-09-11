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
            dict_file = self.get_from_db(ref_type, 'dict_file')
            return {"ref_type": ref_type,
                    "fasta_file": fasta,
                    "dbsnp_file": dbsnp,
                    "dict_file": dict_file}
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
        for idx, f in enumerate([fasta0, fasta1]):
            if f:
                logger.info("Detected fasta: {}".format(f))
            else:
                logger.info("Did not detect fasta {}: {}".format(idx, f))
        return {'fasta0': fasta0, 'fasta1': fasta1}

    ###################################
    # UPLOAD METHODS
    def set_fastas(self, fasta0, fasta1):
        for f in [fasta0, fasta1]:
            continue
            if not os.path.isfile(fasta0):
                raise Exception("Not a valid fasta: {}".format(f))
        self.outputs_dict[self.dataset_name] = {}
        self.outputs_dict[self.dataset_name]['fasta0'] = fasta0
        self.outputs_dict[self.dataset_name]['fasta1'] = fasta1

    def upload_to_db(self, key_name, value):
        """ this dataset must already exist """
        url = "http://data.edicogenome.com/api/set"
        data = {
            'set_x': key_name,
            'name': self.dataset_name,
            'value': value
            }

        logger.info("Uploading to db, Url: {}, Data: {}".format(url, data))
        res = requests.post(url, params=data)
        res.raise_for_status()

    def post_reads(self, fq1_path, fq2_path):
        self.upload_to_db('fastq_location_1', fq1_path)
        self.upload_to_db('fastq_location_2', fq2_path)

    def post_truth_vcf(self, truth_vcf):
        self.upload_to_db('truth_set_vcf', truth_vcf)
        self.upload_to_db('gold_roc_flag', 1)
        # self.upload_to_db('gatk_vcf_pre_filtered', truth_vcf)

    def dataset_create_or_update(self, dataset_name, reference_type):

        logger.info("Check if dataset exists")
        url = "http://data.edicogenome.com/api/get"

        _filter = {'name': dataset_name, 'get_x': 'ref_type'}

        r = requests.get(url, params=_filter)
        r.raise_for_status()
        r = r.text

        if r:
            logger.info("Dataset exists: update existing")
            self.upload_to_db("ref_type", reference_type)
        else:
            logger.info("Dataset does not exist: create new")
            self.create_dna_dataset(dataset_name, reference_type)


    def create_dna_dataset(self, dataset_name, reference_type):
        data = {
            "name": dataset_name,
            "group_id": "2",
            "user": "simulator",
            "values": [
                {"key": 16, "value": 1, "type": "checkbox"}, #  gold set true
                {"key": 37, "value": "", "type": "text"},
                {"key": 38, "value": "", "type": "text"},
                {"key": 39, "value": "", "type": "text"},
                {"key": 40, "value": "", "type": "text"},
                {"key": 41, "value": "", "type": "text"},
                {"key": 42, "value": "", "type": "text"},
                {"key": 43, "value": reference_type, "type": "fk"},
                {"key": 55, "value": "33", "type": "text"},
                {"key": 294, "value": "DRAGEN_RGID", "type": "text"},
                {"key": 295, "value": "DRAGEN_RGSM", "type": "text"},
                {"key": 296, "value": "ILLUMINA", "type": "text"},
                {"key": 297, "value": "1", "type": "text"},
            ],
        }

        headers = {'content-type': 'application/json'}

        r = requests.post("http://data.edicogenome.com/api/dataset/submit",
            data=json.dumps(data), headers=headers)

        try:
            r.raise_for_status()
            logger.info("Dataset updated")
        except Exception as e:
            logger.info("status code: {}".format(r.status_code))
            logger.info("reason: {}".format(r.reason))
            logger.info(e)
            raise

