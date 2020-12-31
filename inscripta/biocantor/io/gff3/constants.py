from enum import Enum

from inscripta.biocantor.util.enum import HasMemberMixin

# in all GFF3 key-value pairs, we escape equals, semicolon, whitespace, ">" and commas
ENCODING_MAP = {"\\\t": "%09", ";": "%3B", "=": "%3D", "\\\n": "%0A", "\\\r": "%0D", ">": "%3E", " ": "%20"}
ENCODING_MAP_WITH_COMMA = {
    "\\\t": "%09",
    ";": "%3B",
    "=": "%3D",
    "\\\n": "%0A",
    "\\\r": "%0D",
    ",": "%2C",
    ">": "%3E",
    " ": "%20",
}
ENCODING_PATTERN = r"({})".format("|".join(k for k in ENCODING_MAP))
ENCODING_PATTERN_WITH_COMMA = r"({})".format("|".join(k for k in ENCODING_MAP_WITH_COMMA))
GFF_SOURCE = "BioCantor"
NULL_COLUMN = "."
ATTRIBUTE_SEPARATOR = ","


class GFF3Headers(Enum):
    HEADER = "##gff-version 3"
    FASTA_HEADER = "##FASTA"
    SEQUENCE_HEADER = "##sequence-region {symbol} 1 {length}"


class _GFF3ReservedQualifiers(HasMemberMixin):
    """These are special reserved qualifiers that BioCantor *does not currently support*.

    All of these, if present in a qualifiers dictionary, will be simply included as-is.

    See :class:`~biocantor.util.gff3.rows.GFFAttributes` for more information."""

    ALIAS = "Alias"
    TARGET = "Target"
    DBXREF = "Dbxref"
    GAP = "Gap"
    DERIVES_FROM = "Derives_from"
    NOTE = "Note"
    ONTOLOGY_TERM = "Ontology_term"


class BioCantorGFF3ReservedQualifiers(HasMemberMixin):
    """This is the subset of GFF3 reserved qualifiers that BioCantor currently reserves"""

    NAME = "Name"
    PARENT = "Parent"
    ID = "ID"


GFF3ReservedQualifiers = HasMemberMixin(
    "GenBankFeatures",
    [[i.name, i.value] for j in [BioCantorGFF3ReservedQualifiers, _GFF3ReservedQualifiers] for i in j],
)


class BioCantorQualifiers(Enum):
    """These are qualifiers that are added when exporting from BioCantor to GFF3, if they exist on the object."""

    TRANSCRIPT_ID = "transcript_id"
    TRANSCRIPT_NAME = "transcript_name"
    TRANSCRIPT_TYPE = "transcript_biotype"
    PROTEIN_ID = "protein_id"
    GENE_ID = "gene_id"
    GENE_NAME = "gene_name"
    GENE_TYPE = "gene_biotype"
    FEATURE_ID = "feature_id"
    FEATURE_SYMBOL = "feature_name"
    FEATURE_TYPE = "feature_type"


class GFF3FeatureTypes(Enum):
    """These are feature types seen in GFF3 files we are parsing that we currently understand."""

    GENE = "gene"
    TRANSCRIPT = "transcript"
    CDS = "CDS"
    EXON = "exon"
    START_CODON = "start_codon"
    STOP_CODON = "stop_codon"
    INTRON = "intron"
    ENHANCER = "enhancer"


BioCantorFeatureTypes = Enum(
    "BioCantorFeatureTypes",
    [
        (f.name, f.value)
        for f in GFF3FeatureTypes
        if f in [GFF3FeatureTypes.GENE, GFF3FeatureTypes.TRANSCRIPT, GFF3FeatureTypes.CDS, GFF3FeatureTypes.EXON]
    ],
)
BioCantorFeatureTypes.__doc__ = "These are the feature types currently supported by BioCantor when writing to GFF3."
