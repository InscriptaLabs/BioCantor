"""
Special feature arithmetic operations for CDSs, codons and translation.

Container classes wrap locations to model genes, transcripts and generic genomic intervals.
"""

from biocantor.gene.biotype import Biotype  # noqa F401
from biocantor.gene.cds_frame import CDSPhase, CDSFrame  # noqa F401
from biocantor.gene.codon import Codon, TranslationTable  # noqa F401
from biocantor.gene.cds import CDSInterval  # noqa F401
from biocantor.gene.feature import FeatureInterval, FeatureIntervalCollection  # noqa F401
from biocantor.gene.transcript import TranscriptInterval  # noqa F401
from biocantor.gene.collections import (  # noqa F401
    AnnotationCollection,
)
from biocantor.gene.gene import GeneInterval  # noqa F401
from biocantor.gene.variants import VariantInterval, VariantIntervalCollection  # noqa F401
