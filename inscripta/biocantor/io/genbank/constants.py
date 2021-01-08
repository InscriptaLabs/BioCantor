"""
GenBank parsing constants. Records common feature types and their relationships as enumerations.
"""

from enum import Enum, IntEnum

from inscripta.biocantor.util.enum import HasMemberMixin


class GenBankParserType(IntEnum):
    """Currently implemented types of GenBank parsing. These types of files are not necessarily mutually exclusive,
    but the parsing implementations are different.
    """

    SORTED = 1
    LOCUS_TAG = 2


class MetadataFeatures(HasMemberMixin):
    """GenBank metadata features BioCantor understands."""

    SOURCE = "source"


class GeneFeatures(HasMemberMixin):
    """GenBank gene features BioCantor understands."""

    GENE = "gene"


class TranscriptFeatures(HasMemberMixin):
    """GenBank transcript features types BioCantor understands."""

    CODING_TRANSCRIPT = "mRNA"
    NONCODING_TRANSCRIPT = "ncRNA"
    TRANSFER_RNA = "tRNA"
    RIBOSOMAL_RNA = "rRNA"
    MISC_RNA = "misc_RNA"
    TM_RNA = "tmRNA"


class IntervalFeatures(HasMemberMixin):
    """GenBank interval features types BioCantor understands."""

    CDS = "CDS"
    EXON = "exon"


GenBankFeatures = HasMemberMixin(
    "GenBankFeatures",
    [[i.name, i.value] for j in [GeneFeatures, TranscriptFeatures, IntervalFeatures, MetadataFeatures] for i in j],
)


class GenbankFlavor(Enum):
    PROKARYOTIC = 1
    EUKARYOTIC = 2
