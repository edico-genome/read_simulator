""" Pipeline classes """
from pipelines_base import PipelinesBase
from modules.modules import *

class VLRD(PipelinesBase):
    modules = [VLRDVCF, Pirs]

class AltContig(PipelinesBase):
    modules = [AltContigVCF, Pirs]  # AltContigPirsTruthSam]

class CNV_Exomes(PipelinesBase):
    modules = [CNVExomeModFastas, Capsim]

class CNV_WHGs(PipelinesBase):
    modules = [RSVSIM_VCF, VCF2Fastas, Pirs]

class Normal(PipelinesBase):
    modules = [CleanVCF, VCF2Fastas, Pirs]

class Tumor(PipelinesBase):
    modules = [TumorHead, VCF2Fastas, Pirs_Tumor]

class Test(PipelinesBase):
    modules = [TestMod]
