""" Pipeline classes """
from pipelines_base import PipelinesBase
from modules.modules import *

class VLRD(PipelinesBase):
    modules = [VLRDVCF, Pirs]

class AltContig(PipelinesBase):
    modules = [AltContigVCF, Pirs]  # AltContigPirsTruthSam]

class CNVExomes(PipelinesBase):
    modules = [CNVExomeModFastas, Capsim]

class CNV_WHGs(PipelinesBase):
    modules = [CNV_WHG_ModFastas, Pirs]
