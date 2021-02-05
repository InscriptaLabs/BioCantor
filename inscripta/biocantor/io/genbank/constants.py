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


class MetadataFeatures(str, HasMemberMixin):
    """GenBank metadata features BioCantor understands."""

    SOURCE = "source"


class GeneFeatures(str, HasMemberMixin):
    """GenBank gene features BioCantor understands."""

    GENE = "gene"


class FeatureCollectionFeatures(str, HasMemberMixin):
    """GenBank feature collections BioCantor understands."""

    FEATURE_COLLECTION = "misc_feature"


class TranscriptFeatures(str, HasMemberMixin):
    """GenBank transcript features types BioCantor understands."""

    CODING_TRANSCRIPT = "mRNA"
    NONCODING_TRANSCRIPT = "ncRNA"
    TRANSFER_RNA = "tRNA"
    RIBOSOMAL_RNA = "rRNA"
    MISC_RNA = "misc_RNA"
    TM_RNA = "tmRNA"


class GeneIntervalFeatures(str, HasMemberMixin):
    """GenBank interval features types BioCantor understands. These do not match

    :class:`~biocantor.io.gff3.constants.BioCantorFeatureTypes` because GenBank has length limitations
    on feature types.
    """

    CDS = "CDS"
    EXON = "exon"


class FeatureIntervalFeatures(str, HasMemberMixin):
    """GenBank interval features types BioCantor understands. These do not match

    :class:`~biocantor.io.gff3.constants.BioCantorFeatureTypes` because GenBank has length limitations
    on feature types.
    """

    FEATURE_INTERVAL = "feat_interval"


class KnownQualifiers(str, Enum):
    """GenBank qualifiers that have special meaning"""

    GENE = "gene"
    LOCUS_TAG = "locus_tag"
    GENE_ID = "gene_id"
    TRANSCRIPT_ID = "transcript_id"
    PROTEIN_ID = "protein_id"
    PRODUCT = "product"
    GENE_NAME = "gene_name"
    GBKEY = "gbkey"
    DBXREF = "db_xref"
    GENE_SYNONYM = "gene_synonym"
    CODON_START = "codon_start"
    FEATURE_CLASS = "feature_class"


# Feature types that mark genes. Used to separate genes from features.
GENBANK_GENE_FEATURES = {"gene", "mRNA", "ncRNA", "tRNA", "rRNA", "misc_RNA", "tmRNA", "CDS", "exon"}


class GenbankFlavor(Enum):
    """GenBank files are formatted differently at NCBI if the species is prokaryotic or eukaryotic. The main
    difference is the presence of a transcript level feature on eukaryotic genomes."""

    PROKARYOTIC = 1
    EUKARYOTIC = 2
