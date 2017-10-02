""" Pipeline classes """

from pipelines_base import PipelinesBase
from modules.modules import *
import logging

logger = logging.getLogger(__name__)


class VLRD(PipelinesBase):
    modules = [VLRDVCF, Pirs]


class AltContig(PipelinesBase):
    # modules = [AltContigVCF, Pirs, AltContigPirsTruthSam]
    modules = [AltContigVCF, Pirs]


class CNVExomes(PipelinesBase):
    modules = [CNVExomeModFastas, Capsim]
