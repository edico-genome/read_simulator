""" Pipeline classes """
from pipelines_base import PipelinesBase
from modules.modules import *

class VLRD(PipelinesBase):
    modules = [FastaTrimmer, BedTrimmer, VLRD_VCF, VCF2Fastas, Pirs, PirsGoldBam]

class VLRD_mason(PipelinesBase):
    modules = [FastaTrimmer, BedTrimmer, VLRD_VCF, VCF2Fastas, Mason]

class AltContig(PipelinesBase):
    modules = [AltContigVCF, Pirs]

class CNV_Exomes(PipelinesBase):
    modules = [PrepCNVprobes, FastaTrimmer, CNVgdbVCF, VCF2Fastas, Capsim]

class CNV_WHGs_Mixed(PipelinesBase):
    modules = [FastaTrimmer, BedTrimmer, CNVgdbVCF, VCF2Fastas, Pirs]

class CNV_WHGs(PipelinesBase):
    modules = [RSVSIM_VCF, ChrTrimmer, VCF2Fastas, Pirs]

class Mason(PipelinesBase):
    modules = [CleanVCF, ChrTrimmer, Mason]

class Normal(PipelinesBase):
    # modules = [CleanVCF, BedTrimmer, VCF2Fastas, Pirs]
    modules = [FastaTrimmer, BedTrimmer, CleanVCF, VCF2Fastas, Pirs]

class Normal_Exome(PipelinesBase):
    # modules = [CleanVCF, BedTrimmer, VCF2Fastas, Pirs]
    modules = [PrepCNVprobes, FastaTrimmer, CleanVCF, VCF2Fastas, Capsim]


class Tumor(PipelinesBase):
    modules = [ChrTrimmer, VCF2BamsurgeonBed, Bamsurgeon]

class Test(PipelinesBase):
    modules = [TestMod]
