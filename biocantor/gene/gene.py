from functools import reduce
from typing import List, Optional, Dict, Hashable, Iterable, Iterator, Any, Union, Set, TYPE_CHECKING, Type
from uuid import UUID

from biocantor import SequenceType
from biocantor.exc import (
    InvalidAnnotationError,
    DuplicateTranscriptError,
    NoncodingTranscriptError,
    NoSuchAncestorException,
)
from biocantor.gene.feature import FeatureInterval, CDSPhase
from biocantor.gene.transcript import TranscriptInterval, Biotype, CDSInterval
from biocantor.gene.biotype import UNKNOWN_BIOTYPE
from biocantor.gene.interval import AbstractFeatureIntervalCollection, IntervalType, QualifierValue
from biocantor.io.gff3.constants import BioCantorQualifiers, GFF_SOURCE, BioCantorFeatureTypes, NULL_COLUMN
from biocantor.io.gff3.exc import GFF3MissingSequenceNameError
from biocantor.io.gff3.rows import GFFRow, GFFAttributes, GTFRow, GTFAttributes
from biocantor.location import Location
from biocantor.parent import Parent
from biocantor.sequence import Sequence
from biocantor.util.bins import bins
from biocantor.util.hashing import digest_object

if TYPE_CHECKING:
    from biocantor.gene.variants import VariantIntervalCollection, VariantInterval


class GeneInterval(AbstractFeatureIntervalCollection):
    """
    A GeneInterval is a collection of :class:`~biocantor.gene.transcript.TranscriptInterval` for a specific locus.

    This is a traditional gene model. By this, I mean that there is one continuous region that defines the gene.
    This region then contains 1 to N subregions that are transcripts. These transcripts may or may not be coding,
    and there is no requirement that all transcripts have the same type. Each transcript consists of one or more
    intervals, and can exist on either strand. There is no requirement that every transcript exist on the same strand.

    The ``Strand`` of this gene interval is always the *plus* strand.

    This cannot be empty; it must have at least one transcript.

    If a ``primary_transcript`` is not provided, then it is inferred by the hierarchy of longest CDS followed by
    longest isoform.
    """

    interval_type = IntervalType.TRANSCRIPT
    _identifiers = ["gene_id", "gene_symbol", "locus_tag"]

    def __init__(
        self,
        transcripts: List[TranscriptInterval],
        guid: Optional[UUID] = None,
        gene_id: Optional[str] = None,
        gene_symbol: Optional[str] = None,
        gene_type: Optional[Biotype] = None,
        locus_tag: Optional[str] = None,
        qualifiers: Optional[Dict[Hashable, List[QualifierValue]]] = None,
        sequence_name: Optional[str] = None,
        sequence_guid: Optional[UUID] = None,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ):
        if not transcripts:
            raise InvalidAnnotationError("GeneInterval must have transcripts")

        self.transcripts = transcripts
        self.gene_id = gene_id
        self.gene_symbol = gene_symbol
        self.gene_type = gene_type
        self.locus_tag = locus_tag
        self.sequence_name = sequence_name
        self.sequence_guid = sequence_guid
        self._parent_or_seq_chunk_parent = parent_or_seq_chunk_parent
        # qualifiers come in as a List, convert to Set
        self._import_qualifiers_from_list(qualifiers)
        self.primary_transcript = AbstractFeatureIntervalCollection._find_primary_feature(self.transcripts)

        # start/end are assumed to be in genomic coordinates, and then _initialize_location
        # will transform them into chunk-relative coordinates if necessary
        self.start = self.genomic_start = min(tx.start for tx in self.transcripts)
        self.end = self.genomic_end = max(tx.end for tx in self.transcripts)
        self._initialize_location(self.start, self.end, parent_or_seq_chunk_parent)
        self.bin = bins(self.start, self.end, fmt="bed")

        if guid is None:
            self.guid = digest_object(
                self.chunk_relative_location,
                self.gene_id,
                self.gene_symbol,
                self.gene_type,
                self.locus_tag,
                self.sequence_name,
                self.qualifiers,
                self.children_guids,
            )
        else:
            self.guid = guid

        self.guid_map = {}
        for tx in self.transcripts:
            if tx.guid in self.guid_map:
                raise DuplicateTranscriptError(f"Guid {tx.guid} found more than once in this GeneInterval")
            self.guid_map[tx.guid] = tx

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(identifiers={self.identifiers}, "
            f"Intervals:{','.join(str(f) for f in self.transcripts)})"
        )

    def iter_children(self) -> Iterable[TranscriptInterval]:
        yield from self.transcripts

    @property
    def is_coding(self) -> bool:
        """One or more coding isoforms?"""
        return any(tx.is_coding for tx in self.transcripts)

    @property
    def id(self) -> str:
        """Returns the ID of this gene. Provides a shared API across genes/transcripts and features."""
        return self.gene_id

    @property
    def name(self) -> str:
        """Returns the name of this gene. Provides a shared API across genes/transcripts and features."""
        return self.gene_symbol

    @property
    def children_guids(self):
        return {x.guid for x in self.transcripts}

    def to_dict(self, chromosome_relative_coordinates: bool = True) -> Dict[str, Any]:
        """Convert to a dict usable by :class:`~biocantor.io.models.GeneIntervalModel`."""
        return dict(
            transcripts=[tx.to_dict(chromosome_relative_coordinates) for tx in self.transcripts],
            gene_id=self.gene_id,
            gene_symbol=self.gene_symbol,
            gene_type=self.gene_type.name if self.gene_type else None,
            locus_tag=self.locus_tag,
            qualifiers=self._export_qualifiers_to_list(),
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            gene_guid=self.guid,
        )

    @staticmethod
    def from_dict(vals: Dict[str, Any], parent_or_seq_chunk_parent: Optional[Parent] = None) -> "GeneInterval":
        """Build a :class:`GeneInterval` from a dictionary representation"""
        return GeneInterval(
            transcripts=[TranscriptInterval.from_dict(x, parent_or_seq_chunk_parent) for x in vals["transcripts"]],
            gene_id=vals["gene_id"],
            gene_symbol=vals["gene_symbol"],
            gene_type=Biotype[vals["gene_type"]] if vals["gene_type"] else None,
            locus_tag=vals["locus_tag"],
            qualifiers=vals["qualifiers"],
            sequence_name=vals["sequence_name"],
            sequence_guid=vals["sequence_guid"],
            guid=vals["gene_guid"],
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    def get_primary_transcript(self) -> Union[TranscriptInterval, None]:
        """Get the primary transcript, if it exists."""

        return self.primary_transcript

    def get_primary_cds(self) -> Union[CDSInterval, None]:
        """Get the CDS of the primary transcript, if it exists."""
        if self.get_primary_transcript() is not None:
            return self.primary_transcript.cds

    def get_primary_transcript_sequence(self) -> Union[Sequence, None]:
        """Get the sequence of the primary transcript, if it exists."""
        if self.get_primary_transcript() is not None:
            return self.primary_transcript.get_spliced_sequence()

    def get_primary_feature(self) -> Union[TranscriptInterval, None]:
        """Convenience function that provides shared API between features and transcripts"""
        return self.get_primary_transcript()

    def get_primary_feature_sequence(self) -> Union[Sequence, None]:
        """Convenience function that provides shared API between features and transcripts"""
        return self.get_primary_transcript_sequence()

    def get_primary_cds_sequence(self) -> Union[Sequence, None]:
        """Get the sequence of the primary transcript, if it exists."""
        if self.get_primary_transcript() is not None:
            return self.primary_transcript.get_cds_sequence()

    def get_primary_protein(self) -> Union[Sequence, None]:
        """Get the protein sequence of the primary transcript, if it exists."""
        if self.get_primary_cds() is not None:
            return self.primary_transcript.get_protein_sequence()

    def _produce_merged_feature(self, intervals: List[Location]) -> FeatureInterval:
        """Wrapper function used by both :func:`GeneInterval.get_merged_transcript`
        and :func:`GeneInterval.get_merged_cds`.
        """
        merged = reduce(lambda x, y: x.union(y), intervals)
        interval_starts = [x.start for x in merged.blocks]
        interval_ends = [x.end for x in merged.blocks]

        return FeatureInterval(
            interval_starts=interval_starts,
            interval_ends=interval_ends,
            strand=self.chunk_relative_location.strand,
            qualifiers=self._export_qualifiers_to_list(),
            sequence_guid=self.sequence_guid,
            sequence_name=self.sequence_name,
            feature_types=[self.gene_type.name],
            feature_name=self.gene_symbol,
            feature_id=self.gene_id,
            guid=self.guid,
            parent_or_seq_chunk_parent=self.chunk_relative_location.parent,
        )

    def get_merged_feature(self) -> FeatureInterval:
        """Generate a single :class:`~biocantor.gene.feature.FeatureInterval` that merges all child features together.

        This inherently has no translation and so is returned as a generic feature, not a transcript.
        """
        return self.get_merged_transcript()

    def get_merged_transcript(self) -> FeatureInterval:
        """Generate a single :class:`~biocantor.gene.feature.FeatureInterval` that merges all child features together.

        This inherently has no translation and so is returned as a generic feature, not a transcript.
        """
        intervals = []
        for tx in self.transcripts:
            for i in tx.chromosome_location.blocks:
                intervals.append(i)
        return self._produce_merged_feature(intervals)

    def get_merged_cds(self) -> FeatureInterval:
        """Generate a single :class:`~biocantor.gene.feature.FeatureInterval` that merges all CDS intervals."""
        intervals = []
        for tx in self.transcripts:
            if tx.is_coding:
                for i in tx.cds.chromosome_location.blocks:
                    intervals.append(i)
        if not intervals:
            raise NoncodingTranscriptError("No CDS transcripts found on this gene")
        return self._produce_merged_feature(intervals)

    def export_qualifiers(self) -> Dict[Hashable, Set[str]]:
        """Exports qualifiers for GFF3/GenBank export"""
        qualifiers = self.qualifiers.copy()
        for key, val in [
            [BioCantorQualifiers.GENE_ID.value, self.gene_id],
            [BioCantorQualifiers.GENE_NAME.value, self.gene_symbol],
            [BioCantorQualifiers.GENE_TYPE.value, self.gene_type.name if self.gene_type else UNKNOWN_BIOTYPE],
            [BioCantorQualifiers.LOCUS_TAG.value, self.locus_tag],
        ]:
            if not val:
                continue
            if key not in qualifiers:
                qualifiers[key] = set()
            qualifiers[key].add(val)
        return qualifiers

    def query_by_guids(self, id_or_ids: Union[UUID, List[UUID]]) -> Optional["GeneInterval"]:
        """Filter this gene interval object by a list of unique IDs.

        Args:
            id_or_ids: List of GUIDs, or unique IDs. Can also be a single ID.

        Returns:
           :class:`GeneInterval`, or None if there are no matching guids.
        """
        if isinstance(id_or_ids, UUID):
            ids = [id_or_ids]
        else:
            ids = id_or_ids

        txs = [self.guid_map[i] for i in ids if i in self.guid_map]
        if txs:
            return GeneInterval(
                transcripts=txs,
                gene_symbol=self.gene_symbol,
                gene_id=self.gene_id,
                gene_type=self.gene_type,
                locus_tag=self.locus_tag,
                qualifiers=self._export_qualifiers_to_list(),
                sequence_name=self.sequence_name,
                sequence_guid=self.sequence_guid,
                guid=self.guid,
                parent_or_seq_chunk_parent=self.chunk_relative_location.parent,
            )

    def _to_gff_or_gtf(
        self,
        chromosome_relative_coordinates: bool = True,
        raise_on_reserved_attributes: Optional[bool] = True,
        row_type: Union[Type[GFFRow], Type[GTFRow]] = GFFRow,
        attribute_type: Union[Type[GFFAttributes], Type[GTFAttributes]] = GFFAttributes,
    ) -> Iterator[Union[GFFRow, GTFRow]]:
        if not self.sequence_name:
            raise GFF3MissingSequenceNameError("Must have sequence names to export to GFF3.")

        if not chromosome_relative_coordinates and not self.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            raise NoSuchAncestorException(
                "Cannot export GFF in relative coordinates without a sequence_chunk ancestor."
            )

        qualifiers = self.export_qualifiers()

        gene_guid = str(self.guid)

        attributes = attribute_type(
            id=gene_guid,
            qualifiers=qualifiers,
            name=self.gene_symbol,
            parent=None,
            raise_on_reserved_attributes=raise_on_reserved_attributes,
        )

        # gene feature
        if row_type == GFFRow:
            row = row_type(
                self.sequence_name,
                GFF_SOURCE,
                BioCantorFeatureTypes.GENE,
                (self.start if chromosome_relative_coordinates else self.chunk_relative_start) + 1,
                self.end if chromosome_relative_coordinates else self.chunk_relative_end,
                NULL_COLUMN,
                self.chunk_relative_location.strand,
                CDSPhase.NONE,
                attributes,
            )
            yield row

        for tx in self.transcripts:
            if row_type == GFFRow:
                yield from tx.to_gff(
                    gene_guid,
                    qualifiers,
                    chromosome_relative_coordinates=chromosome_relative_coordinates,
                    raise_on_reserved_attributes=raise_on_reserved_attributes,
                )
            else:
                yield from tx.to_gtf(
                    gene_guid,
                    qualifiers,
                    chromosome_relative_coordinates=chromosome_relative_coordinates,
                )

    def to_gff(
        self,
        chromosome_relative_coordinates: bool = True,
        raise_on_reserved_attributes: Optional[bool] = True,
    ) -> Iterator[GFFRow]:
        """Produces iterable of :class:`~biocantor.io.gff3.rows.GFFRow` for this gene and its children.

        Args:
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.
            raise_on_reserved_attributes: If ``True``, then GFF3 reserved attributes such as ``ID`` and ``Name`` present
                in the qualifiers will lead to an exception and not a warning.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
        """
        yield from self._to_gff_or_gtf(
            chromosome_relative_coordinates,
            raise_on_reserved_attributes,
            GFFRow,
            GFFAttributes,
        )

    def to_gtf(
        self,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterator[GTFRow]:
        """Produces iterable of :class:`~biocantor.io.gff3.rows.GTFRow` for this gene and its children.

        Args:
            chromosome_relative_coordinates: Output GTF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.

        Yields:
            :class:`~biocantor.io.gff3.rows.GTFRow`

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
        """
        yield from self._to_gff_or_gtf(
            chromosome_relative_coordinates,
            False,
            GTFRow,
            GTFAttributes,
        )

    def incorporate_variants(self, variants: Union["VariantInterval", "VariantIntervalCollection"]) -> "GeneInterval":
        """
        Incorporate all of the variant(s) for an input VariantInterval or VariantIntervalCollection,
        producing a new GeneInterval with those changes incorporated on every child.
        """
        new_transcripts = [tx.incorporate_variants(variants) for tx in self.transcripts]
        if variants.has_sequence:
            new_parent = variants.parent_with_alternative_sequence
        else:
            new_parent = variants.chunk_relative_location.parent
        return GeneInterval(
            new_transcripts,
            guid=None,  # generate a new Interval GUID based on updated data
            gene_id=self.gene_id,
            gene_symbol=self.gene_symbol,
            gene_type=self.gene_type,
            locus_tag=self.locus_tag,
            qualifiers=self._export_qualifiers_to_list(),
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            parent_or_seq_chunk_parent=new_parent,
        )
