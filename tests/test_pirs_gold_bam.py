import pytest, os
from dSim import Settings, main
from lib.common import create_gold_bam_for_pirs
from pprint import pprint


this_dir_path = os.path.dirname(os.path.realpath(__file__))
outdir = "/mnt/archive/sim_data/dSim/vlrd_chr1_normalNoise_varRate0.005/Pirs"
workdir = "/staging/dSim/vlrd_chr1_normalNoise_varRate0.005/Pirs"
liftoverBasename = "/mnt/archive/sim_data/dSim/" + \
                   "vlrd_chr1_normalNoise_varRate0.005/VCF2Fastas/dsim"
read_info = "/mnt/archive/sim_data/dSim/vlrd_chr1_normalNoise_varRate0.005/Pirs/pirs_100_400.read.info.gz"

class Logger:
    def info(self, msg):
        print msg

    def error(self, msg):
        print msg

class SimpleModule:
    def __init__(self):
        self.module_settings = {}
        self.module_settings["outdir"] = outdir
        self.module_settings["workdir"] = workdir
        self.module_settings["liftoverBasename"] = liftoverBasename
        self.module_settings["read_info"] = read_info
        self.module_settings["reference"] = "hs37d5"
        self.logger = Logger()
        

def test_basic_end2end():
    _sm = SimpleModule()
    create_gold_bam_for_pirs(_sm)
    pprint(_sm.module_settings)
