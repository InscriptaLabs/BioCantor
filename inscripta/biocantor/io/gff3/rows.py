"""
Contains information on how to manage GFF row data. Enforces GFF3 specification rules.
"""
import re
from warnings import warn
from typing import Union, Optional, Any, Hashable, Set, Dict
from dataclasses import dataclass

from inscripta.biocantor.io.gff3.constants import (
    ENCODING_MAP,
    ENCODING_PATTERN,
    ENCODING_MAP_WITH_COMMA,
    ENCODING_PATTERN_WITH_COMMA,
    ATTRIBUTE_SEPARATOR,
    BioCantorGFF3ReservedQualifiers,
    GFF3ReservedQualifiers,
    BioCantorFeatureTypes,
)
from inscripta.biocantor.location import Strand
from inscripta.biocantor.gene.cds_frame import CDSPhase
from inscripta.biocantor.io.gff3.exc import GFF3ExportException, ReservedKeyWarning


class GFFAttributes:
    """
    Stores the attributes (column 9) of a GFF row. These attributes are an arbitrary key-value store, but a few rules
    must be enforced. See below for the documentation from the GFF3 spec.

    Of the pre-defined meanings below, currently we only explicitly manage (and require) the following:

    ``ID``
    ``Name``
    ``Parent``

    ``Parent`` is optional, because top level rows do not have a parent. ``Name`` is also optional, because it is not
    a necessary GFF3 item -- it is reserved for human readable names.

    The remaining names are allowed to exist in the qualifiers, but the user will be warned of their presence.
    Since this is case sensitive, this means a warning will happen for ``Dbxref`` but not for ``dbxref``.

    Column 9 tags have predefined meanings:

    * ``ID``: Indicates the unique identifier of the feature. IDs must be unique within the scope of the GFF file.
    * ``Name``: Display name for the feature. This is the name to be displayed to the user. Unlike IDs, there is no
    requirement that the Name be unique within the file.
    * ``Alias``: A secondary name for the feature. It is suggested that this tag be used whenever a secondary identifier
    for the feature is needed, such as locus names and accession numbers. Unlike ID, there is no requirement that Alias
    be unique within the file.
    * ``Parent``: Indicates the parent of the feature. A parent ID can be used to group exons into transcripts,
    transcripts into  genes, and so forth. A feature may have multiple parents. Parent can *only* be used to indicate a
    partof relationship.
    * ``Target``: Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide alignment. The format of
    the value is "target_id start end [strand]", where strand is optional and may be "+" or "-". If the target_id
    contains spaces, they must be escaped as hex escape %20.
    * ``Gap``: The alignment of the feature to the target if the two are not collinear (e.g. contain gaps).
    The alignment format is taken from the CIGAR format described in the Exonerate documentation.
    http://cvsweb.sanger.ac.uk/cgi-bin/cvsweb.cgi/exonerate?cvsroot=Ensembl).
    See the GFF3 specification for more information.
    * ``Derives_from``: Used to disambiguate the relationship between one feature and another when the relationship is a
    temporal one rather than a purely structural "part of" one. This is needed for polycistronic genes. See the GFF3
    specification for more information.
    * ``Note``: A free text note.
    * ``Dbxref``: A database cross reference. See the GFF3 specification for more information.
    * ``Ontology_term``: A cross reference to an ontology term. See the GFF3 specification for more information.

    Multiple attributes of the same type are indicated by separating the values with the comma "," character, as in:

    .. code-block::

        Parent=AF2312,AB2812,abc-3


    Note that attribute names are case sensitive. "Parent" is not the same as "parent".

    All attributes that begin with an uppercase letter are reserved for later use. Attributes that begin with a
    lowercase letter can be used freely by applications. You can stash any semi-structured data into the database by
    using one or more unreserved (lowercase) tags.
    """

    def __init__(
        self,
        id: str,
        qualifiers: Dict[Hashable, Set[Hashable]],
        *,
        name: Optional[str] = None,
        parent: Optional[str] = None,
    ):
        self.id = id
        self.name = name
        self.parent = parent
        self.attributes = qualifiers

        for val in self.attributes.values():
            if not isinstance(val, set):
                raise GFF3ExportException("Attributes dictionary must be a dictionary of sets.")

    def __str__(self):
        """
        Builds a string representation. Handles fixing case where applicable.

        This means joining the key-value pairs with a semicolon, joining the key-values themselves with a equals sign,
        and escaping semicolons, equal signs and tabs in the key or value.

        The GFF3 spec allows commas, and only has you escape ">" and whitespace in Name/ID/Parent, but for simplicity
        and for integration with downstream tools we escape all of them equally here.
        """
        attrs_list = [
            [BioCantorGFF3ReservedQualifiers.ID.value, GFFAttributes.escape_value(self.id, escape_comma=True)]
        ]

        if self.parent is not None:
            attrs_list.append(
                [
                    BioCantorGFF3ReservedQualifiers.PARENT.value,
                    GFFAttributes.escape_value(self.parent, escape_comma=True),
                ]
            )

        if self.name is not None:
            attrs_list.append(
                [BioCantorGFF3ReservedQualifiers.NAME.value, GFFAttributes.escape_value(self.name, escape_comma=True)]
            )

        for key, value_set in sorted(self.attributes.items()):
            if not value_set:
                continue
            if BioCantorGFF3ReservedQualifiers.has_value(key):
                raise GFF3ExportException(f"Found {key} in an attributes dictionary; this is reserved for internal use")
            elif GFF3ReservedQualifiers.has_value(key):
                warn(f"Attribute {key} was seen in the qualifiers, which is a reserved GFF3 key.", ReservedKeyWarning)
                escaped_key = GFFAttributes.escape_key(str(key), lower=False)
            else:
                escaped_key = GFFAttributes.escape_key(str(key), lower=True)
            escaped_vals = [GFFAttributes.escape_value(value, escape_comma=False) for value in value_set]
            escaped_val = ATTRIBUTE_SEPARATOR.join(sorted(escaped_vals))
            attrs_list.append([escaped_key, escaped_val])
        return ";".join(["=".join(pair) for pair in attrs_list])

    @staticmethod
    def _escape_str(item: str) -> str:
        return re.sub(ENCODING_PATTERN, lambda m: ENCODING_MAP.get(m.group(0)), item)

    @staticmethod
    def _escape_str_with_comma(item: str) -> str:
        return re.sub(ENCODING_PATTERN_WITH_COMMA, lambda m: ENCODING_MAP_WITH_COMMA.get(m.group(0)), item)

    @staticmethod
    def escape_key(key: str, lower: Optional[bool] = False) -> str:
        """Key must be escaped for ``[=;\t]``"""
        r = GFFAttributes._escape_str(key)
        return r.lower() if lower else r

    @staticmethod
    def escape_value(value: Any, escape_comma: Optional[bool] = False) -> str:
        """
        Value must be escaped for ``[=;\t]``; make sure value is also not empty.

        Commas must be escaped for reserved attributes like ID and Name.
        """
        value_str = str(value)
        if escape_comma:
            return GFFAttributes._escape_str_with_comma(value_str) if len(value_str) > 0 else "nan"
        else:
            return GFFAttributes._escape_str(value_str) if len(value_str) > 0 else "nan"


@dataclass
class GFFRow:
    """
    Stores the contents of a GFF row. From the spec:

     * Column 1: ``seqid``
    The ID of the landmark used to establish the coordinate system for the current feature. IDs may contain any
    characters, but must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not
    contain unescaped whitespace and must not begin with an unescaped ">".
    To escape a character in this, or any of the other GFF3 fields, replace it with the percent sign followed by its
    hexadecimal representation. For example, ">" becomes "%E3". See URL Encoding (or: 'What are those "%20"
    codes in URLs?') for details.

     * Column 2: ``source``
    The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this
    feature. Typically this is the name of a piece of software, such as "Genescan" or a database name, such as
    "Genbank." In effect, the source is used to extend the feature ontology by adding a qualifier to the type
    creating a new composite type that is a subclass of the type in the type column. It is not necessary to
    specify a source. If there is no source, put a "." (a period) in this field.

    * Column 3: ``type``
    The type of the feature (previously called the "method"). This is constrained to be either: (a) a term from the
    "lite" sequence ontology, SOFA; or (b) a SOFA accession number. The latter alternative is distinguished using
    the syntax SO:000000. This field is required.

     * Columns 4 & 5: ``start`` and ``end``
    The start and end of the feature, in 1-based integer coordinates, relative to the landmark given in column 1.
    Start is always less than or equal to end.
    For zero-length features, such as insertion sites, start equals end and the implied site is to the right of the
    indicated base in the direction of the landmark. These fields are required.

     * Column 6: ``score``
    The score of the feature, a floating point number. As in earlier versions of the format, the semantics of the score
    are ill-defined. It is strongly recommended that E-values be used for sequence similarity features, and that
    P-values be used for ab initio gene prediction features. If there is no score, put a "." (a period) in this field.

     * Column 7: ``strand``
    The strand of the feature. + for positive strand (relative to the landmark), - for minus strand, and . for features
    that are not stranded. In addition, ? can be used for features whose strandedness is relevant, but unknown.

     * Column 8: ``phase``
    For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame.
    The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the
    beginning of this feature to reach the first base of the next codon. In other words, a phase of "0" indicates
    that the next codon begins at the first base of the region described by the current line, a phase of "1" indicates
    that the next codon begins at the second base of this region, and a phase of "2" indicates that the codon begins
    at the third base of this region. This is NOT to be confused with the frame, which is simply start modulo 3.
    If there is no phase, put a "." (a period) in this field.
    For forward strand features, phase is counted from the start field. For reverse strand features, phase is counted
    from the end field. The phase is required for all CDS features.

     * Column 9: ``attributes``
    A list of feature attributes in the format tag=value. Multiple tag=value pairs are separated by semicolons.
    URL escaping rules are used for tags or values containing the following characters: ",=;". Spaces are allowed in
    this field, but tabs must be replaced with the %09 URL escape. This field is not required.

    See :class:`GFFAttributes` for further information on the ``attributes`` column.
    """

    seqid: str
    source: str
    type: BioCantorFeatureTypes
    start: int
    end: int
    score: Union[str, float]
    strand: Strand
    phase: CDSPhase
    attributes: GFFAttributes

    def __str__(self) -> str:
        return "\t".join(
            (
                str(x)
                for x in [
                    self.seqid,
                    self.source,
                    self.type.value,
                    self.start,
                    self.end,
                    self.score,
                    self.strand.to_symbol(),
                    self.phase.to_gff(),
                    str(self.attributes),
                ]
            )
        )
