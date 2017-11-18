""" Pipeline classes """
from pipelines_base import PipelinesBase
from modules.modules import *

class VLRD(PipelinesBase):
    modules = [VLRD_VCF_FASTA, Pirs]

class VLRD_mason(PipelinesBase):
    modules = [FastaTrimmer, BedTrimmer, VLRD_VCF, VCF2Fastas, Mason]

class AltContig(PipelinesBase):
    modules = [AltContigVCF, Pirs]  # AltContigPirsTruthSam]

class CNV_Exomes(PipelinesBase):
    modules = [RSVSIM_VCF, ChrTrimmer, VCF2Fastas, Capsim]

class CNV_WHGs(PipelinesBase):
    modules = [RSVSIM_VCF, ChrTrimmer, VCF2Fastas, Pirs]

class Mason(PipelinesBase):
    modules = [CleanVCF, ChrTrimmer, Mason]

class Normal(PipelinesBase):
    modules = [CleanVCF, ChrTrimmer, VCF2Fastas, Pirs]

class Tumor(PipelinesBase):
    modules = [ChrTrimmer, VCF2BamsurgeonBed, Bamsurgeon]

class Test(PipelinesBase):
    modules = [TestMod]
