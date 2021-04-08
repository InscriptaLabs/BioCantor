import itertools
import re
from enum import Enum

from inscripta.biocantor.util.enum import HasMemberMixin

# in all GFF3 key-value pairs, we escape equals, semicolon, whitespace, ">" and commas, as well as %
ENCODING_MAP = {"\t": "%09", ";": "%3B", "=": "%3D", "\n": "%0A", "\r": "%0D", ">": "%3E", " ": "%20", "%": "%25"}
ENCODING_MAP_WITH_COMMA = {
    "\t": "%09",
    ";": "%3B",
    "=": "%3D",
    "\n": "%0A",
    "\r": "%0D",
    ",": "%2C",
    ">": "%3E",
    " ": "%20",
    "%": "%25",
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
    """These are qualifiers that are added when exporting from BioCantor to GFF3, if they exist on the object.

    Note that this enum does not filter arbitrary qualifier types, but rather exists to map attributes of an
    interval object on to the keys in the GFF3 attributes map.
    """

    TRANSCRIPT_ID = "transcript_id"
    TRANSCRIPT_NAME = "transcript_name"
    TRANSCRIPT_TYPE = "transcript_biotype"
    PROTEIN_ID = "protein_id"
    PRODUCT = "product"
    GENE_ID = "gene_id"
    GENE_SYMBOL = "gene_name"
    GENE_NAME = "gene_name"
    GENE_TYPE = "gene_biotype"
    FEATURE_ID = "feature_id"
    FEATURE_NAME = "feature_name"
    FEATURE_SYMBOL = "feature_name"
    FEATURE_COLLECTION_NAME = "feature_collection_name"
    FEATURE_COLLECTION_ID = "feature_collection_id"
    FEATURE_COLLETION_TYPE = "feature_collection_type"
    FEATURE_TYPE = "feature_type"
    LOCUS_TAG = "locus_tag"


# build a regex of all possible values, case insensitive
BIOCANTOR_QUALIFIERS_REGEX = re.compile(
    r"({})".format(
        "|".join(
            set.union(
                *[
                    {k.name.lower(), k.value}
                    for k in itertools.chain(BioCantorQualifiers, BioCantorGFF3ReservedQualifiers)
                ]
            )
        )
    )
)


class GFF3GeneFeatureTypes(HasMemberMixin):
    """These are feature types seen in GFF3 files we are parsing that we currently understand."""

    GENE = "gene"
    TRANSCRIPT = "transcript"
    CDS = "CDS"
    EXON = "exon"
    PSEUDOGENE = "pseudogene"


class BioCantorFeatureTypes(HasMemberMixin):
    """These are the feature types currently supported by BioCantor when writing genes, transcripts,
    and feature collections to GFF3. When exporting features, the type of the feature is used directly.

    TODO: Feature types should be explicitly linked to Sequence Ontology types. Biological region is a catch-all
        term that matches both the INSDC specification as well as SO:0001411.
    """

    GENE = "gene"
    TRANSCRIPT = "transcript"
    CDS = "CDS"
    EXON = "exon"
    FEATURE_COLLECTION = "biological_region"
    FEATURE_INTERVAL = "feature_interval"
    FEATURE_INTERVAL_REGION = "subregion"
