"""
Collection classes. The data model is structured into two general categories,
transcripts and features. Each of those are wrapped into genes and feature collections,
respectively. These are then wrapped up into one :class:`AnnotationIntervalCollection`.

:class:`AnnotationIntervalCollections` are the topmost class and hold all possible annotations
for a given interval, as well as the place to find their sequence information.

It is useful to think of transcripts/genes as *transcriptional units*, which mean these data structures
model *transcribed sequence*. In contrast, features are *non-transcribed*, and are meant to model things
such as promoters or transcription factor binding sites.

Each object is capable of exporting itself to BED and GFF3.
"""
import itertools
from abc import ABC, abstractmethod
from methodtools import lru_cache
from functools import reduce
from typing import List, Iterable, Any, Dict, Set, Hashable, Optional, Union, Iterator
from uuid import UUID

from inscripta.biocantor.exc import (
    ValidationException,
    NoncodingTranscriptError,
    InvalidAnnotationError,
    InvalidQueryError,
    NoSuchAncestorException,
)
from inscripta.biocantor.gene.biotype import Biotype, UNKNOWN_BIOTYPE
from inscripta.biocantor.gene.cds import CDSInterval
from inscripta.biocantor.gene.cds_frame import CDSPhase
from inscripta.biocantor.gene.feature import FeatureInterval
from inscripta.biocantor.gene.interval import AbstractInterval, QualifierValue, IntervalType
from inscripta.biocantor.gene.transcript import TranscriptInterval
from inscripta.biocantor.io.gff3.constants import GFF_SOURCE, NULL_COLUMN, BioCantorQualifiers, BioCantorFeatureTypes
from inscripta.biocantor.io.gff3.exc import GFF3MissingSequenceNameError
from inscripta.biocantor.io.gff3.rows import GFFRow, GFFAttributes
from inscripta.biocantor.location import Location
from inscripta.biocantor.location.location_impl import SingleInterval, EmptyLocation
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent.parent import Parent, SequenceType
from inscripta.biocantor.sequence import Sequence
from inscripta.biocantor.util.bins import bins
from inscripta.biocantor.util.hashing import digest_object


class AbstractFeatureIntervalCollection(AbstractInterval, ABC):
    """
    Abstract class for holding groups of feature intervals. The two implementations of this class
    model Genes or non-transcribed FeatureCollections.

    These are always on the same sequence, but can be on different strands.
    """

    def __iter__(self):
        """Iterate over all children of this collection"""
        yield from self.iter_children()

    @property
    def strand(self) -> Strand:
        return self.chromosome_location.strand

    @abstractmethod
    def iter_children(self) -> Iterable["AbstractInterval"]:
        """Iterate over the children"""

    @abstractmethod
    def children_guids(self) -> Set[UUID]:
        """Get all of the GUIDs for children.

        Returns: A set of UUIDs
        """

    @abstractmethod
    def query_by_guids(self, ids: List[UUID]) -> "AbstractFeatureIntervalCollection":
        """Filter this collection object by a list of unique IDs.

        Args:
            ids: List of GUIDs, or unique IDs.
        """

    def _reset_parent(self, parent: Optional[Parent] = None) -> None:
        """Reset parent of this collection, and all of its children.

        NOTE: This function modifies this collection in-place, and does not return a new copy. This is different
        behavior than the base function, and is this way because all of the children of this collection are also
        recursively modified.

        NOTE: Using this function presents the risk that you will change the sequence of this interval. There are no
        checks that the new parent provides the same sequence basis as the original parent.

        This overrides :meth:`~biocantor.gene.feature.AbstractInterval.reset_parent()`. The original function
        will remain applied on the leaf nodes.
        """
        self._location = self._location.reset_parent(parent)
        for child in self:
            child._reset_parent(parent)

    def _initialize_location(self, start: int, end: int, parent_or_seq_chunk_parent: Optional[Parent] = None):
        """
        Initialize the location for this collection. Assumes that the start/end coordinates are genome-relative,
        and builds a chunk-relative location for this.

        Args:
            start: genome-relative start
            end: genome-relative end
            parent_or_seq_chunk_parent: A parent that could be null, genome relative, or sequence chunk relative.
        """
        self._location = SingleInterval(start, end, Strand.PLUS)
        if parent_or_seq_chunk_parent:
            if parent_or_seq_chunk_parent.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
                super()._liftover_this_location_to_seq_chunk_parent(parent_or_seq_chunk_parent)
            else:
                self._reset_parent(parent_or_seq_chunk_parent)

    def get_reference_sequence(self) -> Sequence:
        """Returns the *plus strand* sequence for this interval"""
        return self.chunk_relative_location.extract_sequence()

    @staticmethod
    def _find_primary_feature(
        intervals: Union[List[TranscriptInterval], List[FeatureInterval]]
    ) -> Optional[Union[TranscriptInterval, FeatureInterval]]:
        """
        Used in object construction to find the primary feature. Shared between :class:`GeneInterval`
        and :class:`FeatureIntervalCollection`.
        """
        # see if we were given a primary feature
        primary_feature = None
        for i, interval in enumerate(intervals):
            if interval.is_primary_feature:
                if primary_feature:
                    raise ValidationException("Multiple primary features/transcripts found")
                primary_feature = intervals[i]
        # if no primary interval was given, then infer by longest CDS then longest interval
        # if this is a feature, then there is no CDS, so set that value to 0
        if primary_feature is None:
            interval_sizes = sorted(
                (
                    [interval.cds_size if hasattr(interval, "cds_size") else 0, len(interval), i]
                    for i, interval in enumerate(intervals)
                ),
                key=lambda x: (-x[0], -x[1]),
            )
            primary_feature = intervals[interval_sizes[0][2]]
        return primary_feature

    def _liftover_this_location_to_seq_chunk_parent(
        self,
        parent_or_seq_chunk_parent: Parent,
    ):
        """Lift over this collection and all of its children"""
        super()._liftover_this_location_to_seq_chunk_parent(parent_or_seq_chunk_parent)
        for child in self:
            child._liftover_this_location_to_seq_chunk_parent(parent_or_seq_chunk_parent)


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
                raise InvalidAnnotationError(f"Guid {tx.guid} found more than once in this GeneInterval")
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
        """Build an :class:`GeneInterval` from a dictionary representation"""
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

    def get_merged_transcript(self) -> FeatureInterval:
        """Generate a single :class:`~biocantor.gene.feature.FeatureInterval` that merges all exons together.

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

    def query_by_guids(self, ids: List[UUID]) -> Optional["GeneInterval"]:
        """Filter this gene interval object by a list of unique IDs.

        Args:
            ids: List of GUIDs, or unique IDs.

        Returns:
           :class:`GeneInterval`, or None if there are no matching guids.
        """
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

    def to_gff(self, chromosome_relative_coordinates: bool = True) -> Iterable[GFFRow]:
        """Produces iterable of :class:`~biocantor.io.gff3.rows.GFFRow` for this gene and its children.

        Args:
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
        """
        if not self.sequence_name:
            raise GFF3MissingSequenceNameError("Must have sequence names to export to GFF3.")

        if not chromosome_relative_coordinates and not self.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            raise NoSuchAncestorException(
                "Cannot export GFF in relative coordinates without a sequence_chunk ancestor."
            )

        qualifiers = self.export_qualifiers()

        gene_guid = str(self.guid)

        attributes = GFFAttributes(id=gene_guid, qualifiers=qualifiers, name=self.gene_symbol, parent=None)
        row = GFFRow(
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
            yield from tx.to_gff(gene_guid, qualifiers, chromosome_relative_coordinates=chromosome_relative_coordinates)


class FeatureIntervalCollection(AbstractFeatureIntervalCollection):
    """A FeatureIntervalCollection is arbitrary container of intervals.

    This can be thought of to be analogous to a :class:`GeneInterval`, but for non-transcribed
    features that are grouped in some fashion. An example is transcription factor
    binding sites for a specific transcription factor.

    The :class:`~biocantor.location.strand.Strand` of this feature interval collection is always the *plus* strand.

    This cannot be empty; it must have at least one feature interval.
    """

    interval_type = IntervalType.FEATURE
    _identifiers = ["feature_collection_id", "feature_collection_name", "locus_tag"]

    def __init__(
        self,
        feature_intervals: List[FeatureInterval],
        feature_collection_name: Optional[str] = None,
        feature_collection_id: Optional[str] = None,
        feature_collection_type: Optional[str] = None,
        locus_tag: Optional[str] = None,
        sequence_name: Optional[str] = None,
        sequence_guid: Optional[UUID] = None,
        guid: Optional[UUID] = None,
        qualifiers: Optional[Dict[Hashable, List[QualifierValue]]] = None,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ):
        if not feature_intervals:
            raise InvalidAnnotationError("Must have at least one feature interval.")

        self.feature_intervals = feature_intervals
        self.feature_collection_name = feature_collection_name
        self.feature_collection_id = feature_collection_id
        self.locus_tag = locus_tag
        self.sequence_name = sequence_name
        self.sequence_guid = sequence_guid
        # qualifiers come in as a List, convert to Set
        self._import_qualifiers_from_list(qualifiers)
        self.feature_collection_type = feature_collection_type
        self.feature_types = set.union(*[x.feature_types for x in feature_intervals])

        if not self.feature_intervals:
            raise InvalidAnnotationError("FeatureCollection must have features")

        # start/end are assumed to be in genomic coordinates, and then _initialize_location
        # will transform them into chunk-relative coordinates if necessary
        self.start = self.genomic_start = min(f.start for f in self.feature_intervals)
        self.end = self.genomic_end = max(f.end for f in self.feature_intervals)
        self._initialize_location(self.start, self.end, parent_or_seq_chunk_parent)
        self.bin = bins(self.start, self.end, fmt="bed")

        self.primary_feature = AbstractFeatureIntervalCollection._find_primary_feature(self.feature_intervals)

        if guid is None:
            self.guid = digest_object(
                self.chunk_relative_location,
                self.feature_collection_name,
                self.feature_collection_id,
                self.feature_collection_type,
                self.feature_types,
                self.locus_tag,
                self.sequence_name,
                self.qualifiers,
                self.children_guids,
            )
        else:
            self.guid = guid

        self.guid_map = {}
        for feat in self.feature_intervals:
            if feat.guid in self.guid_map:
                raise InvalidAnnotationError(f"Guid {feat.guid} found more than once in this FeatureIntervalCollection")
            self.guid_map[feat.guid] = feat

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(identifiers={self.identifiers}, "
            f"Intervals:{','.join(str(f) for f in self.feature_intervals)})"
        )

    def iter_children(self) -> Iterable[FeatureInterval]:
        """Iterate over all intervals in this collection."""
        yield from self.feature_intervals

    @property
    def is_coding(self) -> bool:
        """Never coding."""
        return False

    @property
    def children_guids(self) -> Set[UUID]:
        return {x.guid for x in self.feature_intervals}

    @property
    def id(self) -> str:
        """Returns the ID of this feature collection. Provides a shared API across genes/transcripts and features."""
        return self.feature_collection_id

    @property
    def name(self) -> str:
        """Returns the name of this feature collection. Provides a shared API across genes/transcripts and features."""
        return self.feature_collection_name

    def get_primary_feature(self) -> Union[FeatureInterval, None]:
        """Get the primary feature, if it exists."""

        return self.primary_feature

    def get_primary_feature_sequence(self) -> Union[Sequence, None]:
        """Convenience function that provides shared API between features and transcripts"""
        if self.get_primary_feature() is not None:
            return self.get_primary_feature().get_spliced_sequence()

    def to_dict(self, chromosome_relative_coordinates: bool = True) -> Dict[str, Any]:
        """Convert to a dict usable by :class:`~biocantor.io.models.FeatureIntervalCollectionModel`."""
        return dict(
            feature_intervals=[feat.to_dict(chromosome_relative_coordinates) for feat in self.feature_intervals],
            feature_collection_name=self.feature_collection_name,
            feature_collection_id=self.feature_collection_id,
            feature_collection_type=self.feature_collection_type,
            locus_tag=self.locus_tag,
            qualifiers=self._export_qualifiers_to_list(),
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            feature_collection_guid=self.guid,
        )

    @staticmethod
    def from_dict(
        vals: Dict[str, Any], parent_or_seq_chunk_parent: Optional[Parent] = None
    ) -> "FeatureIntervalCollection":
        """Build an :class:`FeatureIntervalCollection` from a dictionary representation"""
        return FeatureIntervalCollection(
            feature_intervals=[
                FeatureInterval.from_dict(x, parent_or_seq_chunk_parent) for x in vals["feature_intervals"]
            ],
            feature_collection_name=vals["feature_collection_name"],
            feature_collection_id=vals["feature_collection_id"],
            feature_collection_type=vals["feature_collection_type"],
            locus_tag=vals["locus_tag"],
            qualifiers=vals["qualifiers"],
            sequence_name=vals["sequence_name"],
            sequence_guid=vals["sequence_guid"],
            guid=vals["feature_collection_guid"],
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    def export_qualifiers(self) -> Dict[Hashable, Set[str]]:
        """Exports qualifiers for GFF3/GenBank export"""
        qualifiers = self.qualifiers.copy()
        for key, val in [
            [BioCantorQualifiers.FEATURE_COLLECTION_ID.value, self.feature_collection_id],
            [BioCantorQualifiers.FEATURE_COLLECTION_NAME.value, self.feature_collection_name],
            [BioCantorQualifiers.LOCUS_TAG.value, self.locus_tag],
            [BioCantorQualifiers.FEATURE_COLLETION_TYPE.value, self.feature_collection_type],
        ]:
            if not val:
                continue
            if key not in qualifiers:
                qualifiers[key] = set()
            qualifiers[key].add(val)
        if self.feature_types:
            qualifiers[BioCantorQualifiers.FEATURE_TYPE.value] = self.feature_types
        return qualifiers

    def query_by_guids(self, ids: List[UUID]) -> Optional["FeatureIntervalCollection"]:
        """Filter this feature collection object by a list of unique IDs.

        Args:
            ids: List of GUIDs, or unique IDs.

        Returns:
           :class:`FeatureIntervalCollection`, or None if there are no matching GUIDs.
        """
        features = [self.guid_map[i] for i in ids if i in self.guid_map]
        if features:
            return FeatureIntervalCollection(
                feature_intervals=features,
                feature_collection_name=self.feature_collection_name,
                feature_collection_id=self.feature_collection_id,
                feature_collection_type=self.feature_collection_type,
                locus_tag=self.locus_tag,
                qualifiers=self._export_qualifiers_to_list(),
                sequence_name=self.sequence_name,
                sequence_guid=self.sequence_guid,
                guid=self.guid,
                parent_or_seq_chunk_parent=self.chunk_relative_location.parent,
            )

    def to_gff(self, chromosome_relative_coordinates: bool = True) -> Iterable[GFFRow]:
        """Produces iterable of :class:`~biocantor.io.gff3.rows.GFFRow` for this feature collection and its
        children.

        Args:
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
        """
        if not self.sequence_name:
            raise GFF3MissingSequenceNameError("Must have sequence names to export to GFF3.")

        if not chromosome_relative_coordinates and not self.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            raise NoSuchAncestorException(
                "Cannot export GFF in relative coordinates without a sequence_chunk ancestor."
            )

        qualifiers = self.export_qualifiers()

        feat_group_id = str(self.guid)

        attributes = GFFAttributes(
            id=feat_group_id, qualifiers=qualifiers, name=self.feature_collection_name, parent=None
        )

        row = GFFRow(
            self.sequence_name,
            GFF_SOURCE,
            BioCantorFeatureTypes.FEATURE_COLLECTION,
            (self.start if chromosome_relative_coordinates else self.chunk_relative_start) + 1,
            self.end if chromosome_relative_coordinates else self.chunk_relative_end,
            NULL_COLUMN,
            self.chunk_relative_location.strand,
            CDSPhase.NONE,
            attributes,
        )
        yield row

        for feature in self.feature_intervals:
            yield from feature.to_gff(
                feat_group_id, qualifiers, chromosome_relative_coordinates=chromosome_relative_coordinates
            )


class AnnotationCollection(AbstractFeatureIntervalCollection):
    """An AnnotationCollection is a container to contain :class:`GeneInterval` and
    :class:`FeatureIntervalCollection`.

    Encapsulates all possible annotations for a given interval on a specific source.

    If no start/end points are provided, the interval for this collection is the min/max of the data it contains. The
    interval for an AnnotationCollection is always on the plus strand.

    An AnnotationCollection can be empty (both ``feature_collections`` and ``genes`` can be ``None``).

    The object provided to ``parent_or_seq_chunk_parent`` must have a ``chromosome`` sequence-type in its ancestry,
    and there must be associated sequence. This object should look like the object produced by the function
    :meth:`biocantor.io.parser.seq_to_parent()`, and represent a full chromosome sequence. This will be automatically
    instantiated if you use the constructor method in :class:`biocantor.io.parser.ParsedAnnotationRecord`,
    which will import the sequence from a BioPython ``SeqRecord`` object.

    If you are using file parsers, then if the associated file types have sequence information (GenBank or GFF3+FASTA),
    then the sequences will also be automatically included when the :class:`~biocantor.io.parser.ParsedAnnotationRecord`
    is returned.

    *Object Bounds*: If `start` is provided, `end` must be provided, and vice versa. If neither are provided, and a
    `parent_or_seq_chunk_parent` is provided, then the bounds of this collection will be inferred from that object,
    if possible. If not possible, the bounds of the collection will be the bounds of the child objects associated.

    It is possible to instantiate a :class:`AnnotationCollection` with a ``sequence_chunk`` as well. A
    ``sequence_chunk`` is a slice of a chromosomal sequence that allows operations without loading an entire chromosome
    into memory. The easiest way to produce the parental relationship required for this object to operate on
    ``sequence_chunk`` is to instantiate via the constructor :meth:`biocantor.io.parser.seq_chunk_to_parent()`,
    to which you provide the slice of sequence, the chromosomal start/end positions of that slice, and a sequence name,
    and the returned Parent object will be suitable for passing to this class.
    """

    _identifiers = ["name"]

    def __init__(
        self,
        feature_collections: Optional[List[FeatureIntervalCollection]] = None,
        genes: Optional[List[GeneInterval]] = None,
        name: Optional[str] = None,
        id: Optional[str] = None,
        sequence_name: Optional[str] = None,
        sequence_guid: Optional[UUID] = None,
        sequence_path: Optional[str] = None,
        qualifiers: Optional[Dict[Hashable, QualifierValue]] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        completely_within: Optional[bool] = None,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ):

        self.feature_collections = feature_collections if feature_collections else []
        self.genes = genes if genes else []
        self.sequence_name = sequence_name
        self.sequence_guid = sequence_guid
        self.sequence_path = sequence_path
        self._name = name
        self._id = id
        # qualifiers come in as a List, convert to Set
        self._import_qualifiers_from_list(qualifiers)

        # we store the sequence explicitly, because this is how we can retain sequence information
        # for empty collections
        if parent_or_seq_chunk_parent and parent_or_seq_chunk_parent.sequence:
            self.sequence = parent_or_seq_chunk_parent.sequence
        else:
            self.sequence = None

        if start is None and end is not None:
            raise InvalidAnnotationError("If end is provided, start must also be provided.")
        elif end is None and start is not None:
            raise InvalidAnnotationError("If start is provided, end must also be provided.")
        # must both be unset
        elif start is None:
            # build the start/end coordinates of this collection. Start by looking at the provided parent
            # to use the coordinates provided there.
            if parent_or_seq_chunk_parent and parent_or_seq_chunk_parent.has_ancestor_of_type(SequenceType.CHROMOSOME):
                chrom_parent = parent_or_seq_chunk_parent.first_ancestor_of_type(SequenceType.CHROMOSOME)
                if chrom_parent.location:
                    start = chrom_parent.location.start
                    end = chrom_parent.location.end

            # if we have children, and the above did not work, then use the children
            # cannot infer a range for an empty collection
            if start is None and not self.is_empty:
                start = min(f.start for f in self.iter_children())
                end = max(f.end for f in self.iter_children())

        if start is None and end is None:
            # if we still have nothing, we are empty
            self._location = EmptyLocation()
        else:
            self._initialize_location(start, end, parent_or_seq_chunk_parent)
            self.start = start
            self.end = end
            self.bin = bins(self.start, self.end, fmt="bed")

        self.completely_within = completely_within

        self.guid_map = {x.guid: x for x in self.iter_children()}

        self.guid = digest_object(
            self._location, self.name, self.sequence_name, self.qualifiers, self.completely_within, self.children_guids
        )

    def __repr__(self):
        return f"{self.__class__.__name__}({','.join(str(f) for f in self.iter_children())})"

    def __len__(self):
        return len(self.feature_collections) + len(self.genes)

    @property
    def is_empty(self) -> bool:
        """Is this an empty collection?"""
        return len(self) == 0

    @property
    def children_guids(self) -> set:
        return {x.guid for x in self.iter_children()}

    @lru_cache(maxsize=1)
    @property
    def hierarchical_children_guids(self) -> Dict[UUID, Set[UUID]]:
        """Returns children GUIDs in their hierarchical structure."""
        retval = {}
        for child in self.iter_children():
            if child.guid in retval:
                raise InvalidAnnotationError(f"Found multiple interval collections with the same GUID: {child.guid}")
            retval[child.guid] = child.children_guids
        return retval

    @lru_cache(maxsize=1)
    @property
    def interval_guids_to_collections(self) -> Dict[UUID, Union[GeneInterval, FeatureIntervalCollection]]:
        """
        For example, if this collection had a gene with two transcripts with GUID ABC and 123, and the gene
        had GUID XYZ, this would return:

        .. code-block::

            {
              "ABC": GeneInterval(guid=XYZ),
              "123": GeneInterval(guid=XYZ)
            }


        Returns:
            A map of sub-feature GUIDs to their containing elements.
        """
        retval = {}
        for child in self.iter_children():
            for interval in child.iter_children():
                if interval.guid in retval:
                    raise InvalidAnnotationError(f"Found multiple child intervals with the same GUID: {interval.guid}")
                retval[interval.guid] = child
        return retval

    @property
    def id(self) -> str:
        """Returns the ID of this collection. Provides a shared API across genes/transcripts and features."""
        return self._id

    @property
    def name(self) -> str:
        """Returns the name of this collection. Provides a shared API across genes/transcripts and features."""
        return self._name

    def iter_children(self) -> Iterator[Union[GeneInterval, FeatureIntervalCollection]]:
        """Iterate over all intervals in this collection, in sorted order."""
        chain_iter = itertools.chain(self.genes, self.feature_collections)
        sort_iter = sorted(chain_iter, key=lambda x: x.start)
        yield from sort_iter

    def to_dict(self, chromosome_relative_coordinates: bool = True) -> Dict[str, Any]:
        """Convert to a dict usable by :class:`~biocantor.io.models.AnnotationCollectionModel`."""
        return dict(
            genes=[gene.to_dict(chromosome_relative_coordinates) for gene in self.genes],
            feature_collections=[
                feature.to_dict(chromosome_relative_coordinates) for feature in self.feature_collections
            ],
            name=self.name,
            id=self.id,
            qualifiers=self._export_qualifiers_to_list(),
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            sequence_path=self.sequence_path,
            start=self.start if chromosome_relative_coordinates else self.chunk_relative_start,
            end=self.end if chromosome_relative_coordinates else self.chunk_relative_end,
            completely_within=self.completely_within,
        )

    @staticmethod
    def from_dict(vals: Dict[str, Any], parent_or_seq_chunk_parent: Optional[Parent] = None) -> "AnnotationCollection":
        """Build an :class:`FeatureIntervalCollection` from a dictionary representation"""
        return AnnotationCollection(
            genes=[GeneInterval.from_dict(x, parent_or_seq_chunk_parent) for x in vals["genes"]],
            feature_collections=[
                FeatureIntervalCollection.from_dict(x, parent_or_seq_chunk_parent) for x in vals["feature_collections"]
            ],
            name=vals["name"],
            id=vals["id"],
            qualifiers=vals["qualifiers"],
            sequence_name=vals["sequence_name"],
            sequence_guid=vals["sequence_guid"],
            sequence_path=vals["sequence_path"],
            start=vals["start"],
            end=vals["end"],
            completely_within=vals["completely_within"],
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    def _subset_parent(self, start: int, end: int) -> Optional[Parent]:
        """
        Subset the Parent of this collection to a new interval, building a chunk parent.

        Args:
            start: Genome relative start position.
            end: Genome relative end position.

        Returns:
            A parent, or ``None`` if this location has no parent, or if start == end (empty interval).
        """
        if not self.chunk_relative_location.parent:
            return None
        # edge case for a now null interval
        elif start == end:
            return None
        # edge case -- we are not actually subsetting at all
        if start == self.start and end == self.end:
            return self.chunk_relative_location.parent

        chrom_ancestor = self.lift_over_to_first_ancestor_of_type(SequenceType.CHROMOSOME)

        chunk_relative_start = chrom_ancestor.parent_to_relative_pos(start)

        # handle the edge case where the end is the end of the current chunk
        if end == self.end:
            chunk_relative_end = (
                self.lift_over_to_first_ancestor_of_type(SequenceType.CHROMOSOME).parent_to_relative_pos(end - 1) + 1
            )
        else:
            chunk_relative_end = self.lift_over_to_first_ancestor_of_type(
                SequenceType.CHROMOSOME
            ).parent_to_relative_pos(end)

        seq_subset = self.chunk_relative_location.extract_sequence()[chunk_relative_start:chunk_relative_end]

        parent_id = chrom_ancestor.parent.id
        # TODO: FIXME: handle circular imports by doing this import within the function
        from inscripta.biocantor.io.parser import seq_chunk_to_parent

        return seq_chunk_to_parent(
            str(seq_subset),
            parent_id,
            start,
            end,
            self.chromosome_location.strand,
            self.chunk_relative_location.parent.sequence.alphabet,
        )

    def _build_new_collection_from_query(
        self,
        genes_to_keep: List[GeneInterval],
        features_collections_to_keep: List[FeatureIntervalCollection],
        start: Optional[int],
        end: Optional[int],
        completely_within: Optional[bool],
    ) -> "AnnotationCollection":
        """Convenience function that wraps functionality to build new collections"""
        seq_chunk_parent = self._subset_parent(start, end)
        return AnnotationCollection.from_dict(
            dict(
                feature_collections=[x.to_dict() for x in features_collections_to_keep],
                genes=[x.to_dict() for x in genes_to_keep],
                name=self.name,
                id=self.id,
                sequence_name=self.sequence_name,
                sequence_guid=self.sequence_guid,
                sequence_path=self.sequence_path,
                qualifiers=self._export_qualifiers_to_list(),
                start=start,
                end=end,
                completely_within=completely_within,
            ),
            parent_or_seq_chunk_parent=seq_chunk_parent,
        )

    def query_by_position(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        coding_only: Optional[bool] = False,
        completely_within: Optional[bool] = True,
        expand_location_to_children: Optional[bool] = False,
    ) -> "AnnotationCollection":
        """Filter this annotation collection object based on positions, sequence, and boolean flags.

        In all cases, the comparisons are made without considering strand.  Intronic queries are still valid.
        In other words, a query from ``[10,20]`` would still return a transcript whose intervals were
        ``[0,9], [21, 30]``.

        The resulting :class:`AnnotationCollection` returned will have a `._location` member whose bounds
        exactly match the query. If ``expand_location_to_children`` is ``True``, then the
        child genes/feature collections will potentially extend beyond this range,
        in order to encapsulate their full length. The resulting gene/feature collections will potentially have a
        reduced set of transcripts/features, if those transcripts/features are outside the query range.
        However, if ``expand_location_to_children`` is ``False``, then the child genes/feature collections
        will have location objects that represent the exact bounds of the query, which means that they
        may be sliced down. If the sliced down coordinates are entirely intronic for any isoform, then
        this isoform will have an EmptyLocation `chunk_relative_location` member, because it is no longer
        possible to have a relationship to the location object associated with this collection.

        Here is an example (equals are exons, dashes are introns):

        .. code-block::

                          10      15      20      25      30      35      40
            Gene1: Tx1:     12============20
                   Tx2:     12======16-17=20--22==25
            Fc1:    F1:     12====15
                    F2:     12======16-17=20--22==25
                    F3:                                           35======40


        Results:

        +------------+------------+--------------------+------------------+
        | start      | end        | completely_within  | result           |
        +============+============+====================+==================+
        | 21         | 22         | True               | EmptyCollection  |
        +------------+------------+--------------------+------------------+
        | 21         | 22         | False              | Tx1,Tx2,F2       |
        +------------+------------+--------------------+------------------+
        | 28         | 35         | False              | EmptyCollection  |
        +------------+------------+--------------------+------------------+
        | 28         | 36         | False              | F3               |
        +------------+------------+--------------------+------------------+
        | 27         | 36         | False              | Tx1,F3           |
        +------------+------------+--------------------+------------------+
        | 24         | 36         | False              | Tx1,Tx2,F2,F3    |
        +------------+------------+--------------------+------------------+

        Args:
            start: Genome relative start position. If not set, will be 0.
            end: Genome relative end position. If not set, will be unbounded.
            coding_only: Filter for coding genes only?
            completely_within: Strict *query* boundaries? If ``False``, features that partially overlap
                will be included in the output. Bins optimization cannot be used, so these queries are slower.
            expand_location_to_children: Should the underlying location objects be expanded so that no
                child gene/transcripts get sliced? If this is ``False``, then the constituent objects may not
                actually represent their full lengths, although the original position information is retained.

        Returns:
           :class:`AnnotationCollection` that may be empty, and otherwise will contain new copies of every
            constituent member.

        Raises:
            InvalidQueryError: If the start/end bounds are not valid. This could be because they exceed the
            bounds of the current interval. It could also happen if ``expand_location_to_children`` is ``True``
            and the new expanded range would exceed the range of an associated sequence chunk.
        """
        # bins are only valid if we have start, end and completely_within
        if completely_within and start and end:
            my_bins = bins(start, end, fmt="bed", one=False)
        else:
            my_bins = None

        # after bins were decided, we can now force start/end to min/max values
        # for exact checking
        start = self.start if start is None else start
        end = self.end if end is None else end
        if start < 0:
            raise InvalidQueryError("Start must be positive")
        elif start > end:
            raise InvalidQueryError("Start must be less than or equal to end")
        elif start < self.start:
            raise InvalidQueryError(
                f"Start {start} must be within bounds of current interval [{self.start}-{self.end})"
            )
        elif end > self.end:
            raise InvalidQueryError(f"End {end} must be within bounds of current interval [{self.start}-{self.end})")
        elif start == end:
            raise InvalidQueryError("Cannot query a 0bp interval (start must not be the same as end).")

        # coordinate_fn will be applied when filtering specific transcripts/features
        query_loc = SingleInterval(start, end, Strand.PLUS, parent=self.chromosome_location.parent)
        if completely_within:
            coordinate_fn = query_loc.contains
        else:
            coordinate_fn = query_loc.has_overlap

        genes_to_keep = []
        features_collections_to_keep = []
        for gene_or_feature_collection in self.iter_children():

            if coding_only and not gene_or_feature_collection.is_coding:
                continue

            # my_bins only exists if completely_within, start and end
            # if no children match these bins, skip
            elif my_bins and not any(child.bin in my_bins for child in gene_or_feature_collection):
                continue

            # regardless of completely_within flag, first just look for overlaps on the gene/feature collection level
            elif coordinate_fn(gene_or_feature_collection.chromosome_location, match_strand=False, full_span=True):
                if gene_or_feature_collection.interval_type == IntervalType.FEATURE:
                    features_collections_to_keep.append(gene_or_feature_collection)
                else:
                    genes_to_keep.append(gene_or_feature_collection)

        # if completely within is False, expand the range of seq_chunk_parent to retain the full span
        # of all child intervals. This prevents features getting cut in half.
        if expand_location_to_children is True and completely_within is False:
            for g_or_fc in itertools.chain(features_collections_to_keep, genes_to_keep):
                if g_or_fc.start < start:
                    start = g_or_fc.start
                if g_or_fc.end > end:
                    end = g_or_fc.end

        # if there is a sequence chunk, then some validation checks must be performed
        if self.chunk_relative_location.parent and self.chunk_relative_location.parent.sequence:
            # not possible to expand range if it exceeds parent bounds
            if start < self.start or end > self.end:
                raise InvalidQueryError(
                    f"Cannot expand range of location to {start}-{end} because the associated sequence chunk "
                    f"lies from {self.start}-{self.end}"
                )

        return self._build_new_collection_from_query(
            genes_to_keep, features_collections_to_keep, start, end, completely_within
        )

    def _return_collection_for_id_queries(
        self, genes_to_keep: List[GeneInterval], features_collections_to_keep: List[FeatureIntervalCollection]
    ) -> "AnnotationCollection":
        """Convenience function shared by all functions that query by identifiers or GUIDs."""

        if genes_to_keep or features_collections_to_keep:
            start = min(self.start, min(x.start for x in itertools.chain(genes_to_keep, features_collections_to_keep)))
            end = max(self.end, max(x.end for x in itertools.chain(genes_to_keep, features_collections_to_keep)))
        else:
            start = self.start
            end = self.end

        return self._build_new_collection_from_query(
            genes_to_keep, features_collections_to_keep, start, end, self.completely_within
        )

    def query_by_guids(self, ids: List[UUID]) -> "AnnotationCollection":
        """Filter this annotation collection object by a list of unique IDs.

        Args:
            ids: List of GUIDs, or unique IDs.

        NOTE: If the children of this collection have GUID collisions, either across genes or features or
        within genes and features, this function will return all members with the matching GUID.

        Returns:
           :class:`AnnotationCollection` that may be empty.
        """
        genes_to_keep = []
        features_collections_to_keep = []
        for i in ids:
            gene_or_feature_collection = self.guid_map.get(i)
            if not gene_or_feature_collection:
                continue
            elif gene_or_feature_collection.interval_type == IntervalType.FEATURE:
                features_collections_to_keep.append(gene_or_feature_collection)
            else:
                genes_to_keep.append(gene_or_feature_collection)

        return self._return_collection_for_id_queries(genes_to_keep, features_collections_to_keep)

    def query_by_interval_guids(self, id_or_ids: Union[UUID, List[UUID]]) -> "AnnotationCollection":
        """Filter this annotation collection object by a list of unique *interval* IDs.

        This function wraps the ``query_by_guid`` function of child GeneInterval/FeatureIntervalCollection
        objects.

        NOTE: If the children of this collection have GUID collisions, either across genes or features or
        within genes and features, this function will return all members with the matching GUID.

        Args:
            id_or_ids: List of GUIDs, or unique IDs. Can also be a single ID.

        Returns:
           :class:`AnnotationCollection` that may be empty.
        """
        if isinstance(id_or_ids, UUID):
            ids = [id_or_ids]
        else:
            ids = id_or_ids

        genes_to_keep = []
        features_collections_to_keep = []
        for child in self.iter_children():
            gene_or_feature_collection = child.query_by_guids(ids)
            if not gene_or_feature_collection:
                continue
            elif gene_or_feature_collection.interval_type == IntervalType.FEATURE:
                features_collections_to_keep.append(gene_or_feature_collection)
            else:
                genes_to_keep.append(gene_or_feature_collection)

        return self._return_collection_for_id_queries(genes_to_keep, features_collections_to_keep)

    def query_by_transcript_interval_guids(self, id_or_ids: Union[UUID, List[UUID]]) -> "AnnotationCollection":
        """Filter this annotation collection object by a list of unique *TranscriptInterval* IDs.

        This function wraps the ``query_by_guid`` function of child GeneInterval objects.

        NOTE: If the children of this collection have GUID collisions, either across genes or features or
        within genes and features, this function will return all members with the matching GUID.

        Args:
            id_or_ids: List of GUIDs, or unique IDs. Can also be a single ID.

        Returns:
           :class:`AnnotationCollection` that may be empty.
        """
        if isinstance(id_or_ids, UUID):
            ids = [id_or_ids]
        else:
            ids = id_or_ids

        genes_to_keep = []
        for child in self.iter_children():
            gene_or_feature_collection = child.query_by_guids(ids)
            if not gene_or_feature_collection:
                continue
            elif gene_or_feature_collection.interval_type == IntervalType.TRANSCRIPT:
                genes_to_keep.append(gene_or_feature_collection)

        return self._return_collection_for_id_queries(genes_to_keep, [])

    def query_by_feature_interval_guids(self, id_or_ids: Union[UUID, List[UUID]]) -> "AnnotationCollection":
        """Filter this annotation collection object by a list of unique *interval* IDs.

        This function wraps the ``query_by_guid`` function of child FeatureIntervalCollection objects.

        NOTE: If the children of this collection have GUID collisions, either across genes or features or
        within genes and features, this function will return all members with the matching GUID.

        Args:
            id_or_ids: List of GUIDs, or unique IDs. Can also be a single ID.

        Returns:
           :class:`AnnotationCollection` that may be empty.
        """
        if isinstance(id_or_ids, UUID):
            ids = [id_or_ids]
        else:
            ids = id_or_ids

        features_collections_to_keep = []
        for child in self.iter_children():
            gene_or_feature_collection = child.query_by_guids(ids)
            if not gene_or_feature_collection:
                continue
            elif gene_or_feature_collection.interval_type == IntervalType.FEATURE:
                features_collections_to_keep.append(gene_or_feature_collection)

        return self._return_collection_for_id_queries([], features_collections_to_keep)

    def query_by_feature_identifiers(self, id_or_ids: Union[str, List[str]]) -> "AnnotationCollection":
        """Filter this annotation collection object by a list of identifiers, or a single identifier.

        Identifiers are not necessarily unique; if your identifier matches more than one interval,
        all matching intervals will be returned. These ambiguous results will be adjacent in the resulting collection,
        but are not grouped or signified in any way.

        This method is ``O(n_ids * m_identifiers)``.

        Args:
            id_or_ids: List of identifiers, or a single identifier.

        Returns:
           :class:`AnnotationCollection` that may be empty.
        """
        if isinstance(id_or_ids, str):
            ids = {id_or_ids}
        else:
            ids = set(id_or_ids)

        genes_to_keep = []
        features_collections_to_keep = []
        for gene_or_feature in self.iter_children():
            if ids & gene_or_feature.identifiers:
                if isinstance(gene_or_feature, FeatureIntervalCollection):
                    features_collections_to_keep.append(gene_or_feature)
                else:
                    genes_to_keep.append(gene_or_feature)

        return self._return_collection_for_id_queries(genes_to_keep, features_collections_to_keep)

    def get_children_by_type(self, child_type: str) -> Union[List[GeneInterval], List[FeatureIntervalCollection]]:
        if child_type.lower() == IntervalType.FEATURE:
            return self.feature_collections
        elif child_type.lower() == IntervalType.TRANSCRIPT:
            return self.genes
        else:
            raise InvalidQueryError(f"Cannot get children of type {child_type}")

    def _unsorted_gff_iter(self, chromosome_relative_coordinates: bool = True) -> Iterable[GFFRow]:
        """Produces iterable of :class:`~biocantor.io.gff3.rows.GFFRow` for this annotation collection and its
        children.

        The positions of the genes will be ordered by genomic position, but may not be globally position sorted
        because it could be the case that children gene/features will overlap. This private function
        exists to provide an iterator to sort in the main ``to_gff()`` function.

        Args:
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`
        """
        for item in self.iter_children():
            yield from item.to_gff(chromosome_relative_coordinates=chromosome_relative_coordinates)

    def to_gff(self, chromosome_relative_coordinates: bool = True) -> Iterable[GFFRow]:
        """Produces iterable of :class:`~biocantor.io.gff3.rows.GFFRow` for this annotation collection and its
        children.

        Args:
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
        """
        yield from sorted(self._unsorted_gff_iter(chromosome_relative_coordinates), key=lambda x: x.start)
