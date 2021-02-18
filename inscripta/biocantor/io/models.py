"""
Data models. These models allow for validation of inputs to a BioCantor model, acting as a JSON schema for serializing
and deserializing the models.
"""
from typing import List, Optional, ClassVar, Type, Dict, Union
from uuid import UUID

from inscripta.biocantor.exc import (
    InvalidCDSIntervalError,
    ValidationException,
    NullSequenceException,
)
from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds import CDSFrame, CDSInterval
from inscripta.biocantor.gene.collections import GeneInterval, FeatureIntervalCollection, AnnotationCollection
from inscripta.biocantor.gene.feature import FeatureInterval
from inscripta.biocantor.gene.transcript import TranscriptInterval
from inscripta.biocantor.location.location import Location
from inscripta.biocantor.location.location_impl import CompoundInterval, SingleInterval
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import Parent
from marshmallow import Schema  # noqa: F401
from marshmallow_dataclass import dataclass


@dataclass
class BaseModel:
    """Base for all of the models."""

    Schema: ClassVar[Type[Schema]] = Schema  # noqa: F811

    class Meta:
        ordered = True

    @staticmethod
    def initialize_location(
        starts: List[int],
        ends: List[int],
        strand: Strand,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ) -> Location:
        """
        Initialize the :class:`Location` object for this interval.

        Args:
            starts: Start positions relative to the chromosome.
            ends: End positions relative to the chromosome.
            strand: Strand relative to the chromosome.
            parent_or_seq_chunk_parent: An optional parent, either as a full chromosome or as a sequence chunk.
        """
        if len(starts) != len(ends):
            raise ValidationException("Number of interval starts does not match number of interval ends")
        elif len(starts) == len(ends) == 1:
            location = SingleInterval(starts[0], ends[0], strand)
        else:
            location = CompoundInterval(starts, ends, strand)
        return BaseModel.liftover_location_to_seq_chunk_parent(location, parent_or_seq_chunk_parent)

    @staticmethod
    def liftover_location_to_seq_chunk_parent(
        location: Location,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ) -> Location:
        """
        BioCantor supports constructing any of the interval classes from a subset of the chromosome. In order to
        be able to set up the coordinate relationship and successfully pull down sequence, this function
        lifts the coordinates from the original annotation object on to this new coordinate system.

        .. code:: python

            parent_1_15 = Parent(
                sequence=Sequence(
                    genome2[1:15],
                    Alphabet.NT_EXTENDED_GAPPED,
                    type="sequence_chunk",
                    parent=Parent(
                        location=SingleInterval(1, 15, Strand.PLUS,
                                               parent=Parent(id="genome_1_15", sequence_type="chromosome"))
                    ),
                )
            )

        Alternatively, if the sequence is coming straight from a file, it will be a :class:`Parent` with a
        :class:`Sequence` attached:

        .. code:: python

            parent = Parent(id="chr1", sequence=Sequence(genome, Alphabet.NT_STRICT, type="chromosome"))

        This convenience function detects which kind of parent is given, and sets up the appropriate location.

        Args:
            location: A location object, likely produced by :meth:`initialize_location()`.
            parent_or_seq_chunk_parent: An optional parent, either as a full chromosome or as a sequence chunk.

        Returns:
            A :class:`Location` object.

        Raises:
            ValidationException: If ``parent_or_seq_chunk_parent`` has no ancestor of type ``chromosome`` or
                ``sequence_chunk``.
            NullSequenceException: If ``parent_or_seq_chunk_parent`` has no usable sequence ancestor.
        """
        if parent_or_seq_chunk_parent is None:
            return location

        elif parent_or_seq_chunk_parent.has_ancestor_of_type("sequence_chunk"):

            chunk_parent = parent_or_seq_chunk_parent.first_ancestor_of_type("sequence_chunk")
            if not chunk_parent.sequence:
                raise NullSequenceException("Must have a sequence if a parent is provided.")

            location = location.reset_parent(chunk_parent.parent)
            sequence_chunk = chunk_parent.sequence
            interval_location_rel_to_chunk = sequence_chunk.location_on_parent.parent_to_relative_location(location)
            interval_rel_to_chunk = interval_location_rel_to_chunk.reset_parent(parent_or_seq_chunk_parent)

            if len(interval_rel_to_chunk) != len(location):
                raise ValidationException(
                    f"Lifted location is no longer the same length as the input "
                    f"(({len(interval_rel_to_chunk)} != ({len(location)})."
                )

            return interval_rel_to_chunk

        # since this is a whole genome, we don't need to lift anything up
        elif parent_or_seq_chunk_parent.has_ancestor_of_type("chromosome"):
            if not parent_or_seq_chunk_parent.first_ancestor_of_type("chromosome").sequence:
                raise NullSequenceException("Must have a sequence if a parent is provided.")

            return location.reset_parent(parent_or_seq_chunk_parent)

        else:
            raise ValidationException("Provided Parent has no sequence of type 'chromosome' or 'sequence_chunk'.")


@dataclass
class FeatureIntervalModel(BaseModel):
    """Data model that allows construction of a :class:`~biocantor.gene.feature.FeatureInterval` object.

    FeatureIntervals can have more than one type, and these types are arbitrary and not controlled by a biotype
    ontology.
    """

    interval_starts: List[int]
    interval_ends: List[int]
    strand: Strand
    qualifiers: Optional[Dict[str, List[Union[str, int, bool, float]]]] = None
    sequence_name: Optional[str] = None
    sequence_guid: Optional[UUID] = None
    feature_interval_guid: Optional[UUID] = None
    feature_guid: Optional[UUID] = None
    feature_types: Optional[List[str]] = None
    feature_name: Optional[str] = None
    feature_id: Optional[str] = None
    is_primary_feature: Optional[bool] = None

    def to_feature_interval(
        self,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ) -> "FeatureInterval":
        """Construct a :class:`~biocantor.gene.feature.FeatureInterval` from a :class:`FeatureIntervalModel`.

        A :class:`~biocantor.parent.Parent` can be provided to allow the sequence-retrieval functions to work.
        """

        location = BaseModel.initialize_location(
            self.interval_starts,
            self.interval_ends,
            self.strand,
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

        return FeatureInterval(
            location=location,
            qualifiers=self.qualifiers,
            sequence_guid=self.sequence_guid,
            sequence_name=self.sequence_name,
            feature_types=self.feature_types,
            feature_name=self.feature_name,
            feature_id=self.feature_id,
            guid=self.feature_interval_guid,
            feature_guid=self.feature_guid,
            is_primary_feature=self.is_primary_feature,
        )

    @staticmethod
    def from_feature_interval(feature_interval: FeatureInterval) -> "FeatureIntervalModel":
        """Convert a :class:`~biocantor.gene.feature.FeatureInterval` to a :class:`FeatureIntervalModel`"""
        return FeatureIntervalModel.Schema().load(feature_interval.to_dict())


@dataclass
class TranscriptIntervalModel(BaseModel):
    """
    Data model that allows construction of a :class:`~biocantor.gene.feature.TranscriptInterval` object.
    """

    exon_starts: List[int]
    exon_ends: List[int]
    strand: Strand
    cds_starts: Optional[List[int]] = None
    cds_ends: Optional[List[int]] = None
    cds_frames: Optional[List[CDSFrame]] = None
    qualifiers: Optional[Dict[str, List[Union[str, int, bool, float]]]] = None
    is_primary_tx: Optional[bool] = None
    transcript_id: Optional[str] = None
    protein_id: Optional[str] = None
    transcript_symbol: Optional[str] = None
    transcript_type: Optional[Biotype] = None
    sequence_name: Optional[str] = None
    sequence_guid: Optional[UUID] = None
    transcript_interval_guid: Optional[UUID] = None
    transcript_guid: Optional[UUID] = None

    def to_transcript_interval(
        self,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ) -> "TranscriptInterval":
        """Construct a :class:`~biocantor.gene.transcript.TranscriptInterval` from a :class:`TranscriptIntervalModel`.

        A :class:`~biocantor.parent.Parent can be provided to allow the sequence-retrieval functions to work.
        """

        if self.cds_starts is not None and self.cds_ends is None:
            raise InvalidCDSIntervalError("If CDS start is defined, CDS end must be defined")
        elif self.cds_starts is None and self.cds_ends is not None:
            raise InvalidCDSIntervalError("If CDS end is defined, CDS start must be defined")
        elif self.cds_starts is not None and self.cds_ends is not None:  # must be coding
            if len(self.cds_starts) != len(self.cds_ends):
                raise InvalidCDSIntervalError("Number of CDS starts does not number of CDS ends")
            elif self.cds_starts[0] < self.exon_starts[0]:
                raise InvalidCDSIntervalError("CDS start must be greater than or equal to exon start")
            elif self.cds_ends[-1] > self.exon_ends[-1]:
                raise InvalidCDSIntervalError("CDS end must be less than or equal to than exon end")
            elif self.cds_frames is None:
                raise InvalidCDSIntervalError("If CDS interval is defined, CDS frames must be defined")
            elif len(self.cds_frames) != len(self.cds_starts):
                raise InvalidCDSIntervalError("Number of CDS frames must match number of CDS starts/ends")

            cds_interval = BaseModel.initialize_location(
                self.cds_starts, self.cds_ends, self.strand, parent_or_seq_chunk_parent=parent_or_seq_chunk_parent
            )
            if not cds_interval:
                raise InvalidCDSIntervalError("CDS must have a non-zero length")

            # length validation happens in the CDS constructor
            cds = CDSInterval(cds_interval, self.cds_frames)
        else:
            cds = None

        location = BaseModel.initialize_location(
            self.exon_starts,
            self.exon_ends,
            self.strand,
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

        return TranscriptInterval(
            location=location,
            cds=cds,
            guid=self.transcript_interval_guid,
            transcript_guid=self.transcript_guid,
            qualifiers=self.qualifiers,
            is_primary_tx=self.is_primary_tx,
            transcript_id=self.transcript_id,
            transcript_symbol=self.transcript_symbol,
            transcript_type=self.transcript_type if self.transcript_type else None,
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            protein_id=self.protein_id,
        )

    @staticmethod
    def from_transcript_interval(transcript_interval: TranscriptInterval) -> "TranscriptIntervalModel":
        """Convert to a :class:`~biocantor.models.TranscriptIntervalModel`"""

        return TranscriptIntervalModel.Schema().load(transcript_interval.to_dict())


@dataclass
class GeneIntervalModel(BaseModel):
    """
    Data model that allows construction of :class:`~biocantor.gene.collections.GeneInterval` object.

    This is a container for one or more :class:`~biocantor.gene.transcript.TranscriptInterval` objects.

    Has additional keys to help query the existing Gene table to see if this gene is there already.
    """

    transcripts: List[TranscriptIntervalModel]
    gene_id: Optional[str] = None
    gene_symbol: Optional[str] = None
    gene_type: Optional[Biotype] = None
    locus_tag: Optional[str] = None
    qualifiers: Optional[Dict[str, List[Union[str, int, bool, float]]]] = None
    sequence_name: Optional[str] = None
    sequence_guid: Optional[UUID] = None
    gene_guid: Optional[UUID] = None

    def to_gene_interval(self, parent_or_seq_chunk_parent: Optional[Parent] = None) -> GeneInterval:
        """Produce a :class:`~biocantor.gene.collections.GeneInterval` from a :class:`GeneIntervalModel`.

        This is the primary method of constructing a :class:`biocantor.gene.collections.GeneInterval`.
        """

        transcripts = [tx.to_transcript_interval(parent_or_seq_chunk_parent) for tx in self.transcripts]

        return GeneInterval(
            transcripts=transcripts,
            guid=self.gene_guid,
            gene_id=self.gene_id,
            gene_symbol=self.gene_symbol,
            gene_type=self.gene_type,
            locus_tag=self.locus_tag,
            qualifiers=self.qualifiers,
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    @staticmethod
    def from_gene_interval(gene: GeneInterval) -> "GeneIntervalModel":
        return GeneIntervalModel.Schema().load(gene.to_dict())


@dataclass
class FeatureIntervalCollectionModel(BaseModel):
    """
    Data model that allows construction of :class:`~biocantor.gene.collections.FeatureCollection` object.

    This is a container for one or more :class:`~biocantor.gene.feature.FeatureInterval` objects.

    Feature Collections do not have a type, but rather their type is the union of all of their child types.
    """

    feature_intervals: List[FeatureIntervalModel]
    feature_collection_name: Optional[str] = None
    feature_collection_id: Optional[str] = None
    locus_tag: Optional[str] = None
    feature_collection_type: Optional[str] = None
    sequence_name: Optional[str] = None
    sequence_guid: Optional[UUID] = None
    feature_collection_guid: Optional[UUID] = None
    qualifiers: Optional[Dict[str, List[Union[str, int, bool, float]]]] = None

    def to_feature_collection(self, parent_or_seq_chunk_parent: Optional[Parent] = None) -> FeatureIntervalCollection:
        """Produce a feature collection from a :class:`FeatureIntervalCollectionModel`."""

        feature_intervals = [feat.to_feature_interval(parent_or_seq_chunk_parent) for feat in self.feature_intervals]

        return FeatureIntervalCollection(
            feature_intervals=feature_intervals,
            feature_collection_name=self.feature_collection_name,
            feature_collection_id=self.feature_collection_id,
            locus_tag=self.locus_tag,
            feature_collection_type=self.feature_collection_type,
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            qualifiers=self.qualifiers,
            guid=self.feature_collection_guid,
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    @staticmethod
    def from_feature_collection(feature_collection: FeatureIntervalCollection) -> "FeatureIntervalCollectionModel":
        return FeatureIntervalCollectionModel.Schema().load(feature_collection.to_dict())


@dataclass
class AnnotationCollectionModel(BaseModel):
    """
    Data model that allows construction of :class:`~biocantor.gene.collections.AnnotationCollection` object.

    This is the highest level type of container, and is most often used to store everything for a genomic
    interval query.

    This container has optional ``start`` and ``end`` members because it can be the product of genomic interval queries
    that are larger than the contents of the object. If those values are not set, when an
    :class:~biocantor.gene.collections.AnnotationCollection` is instantiated, they will be inferred.

    Additionally, this container has an optional ``completely_within`` member that determines if a range query that
    produced this container was done using the ``completely_within`` flag. If this is ``True``, then it may be the case
    that ``start`` is larger than the smallest start position of a member of this collection, and vice versa for
    ``end``.
    """

    feature_collections: Optional[List[FeatureIntervalCollectionModel]] = None
    genes: Optional[List[GeneIntervalModel]] = None
    name: Optional[str] = None
    id: Optional[str] = None
    sequence_name: Optional[str] = None
    sequence_guid: Optional[UUID] = None
    sequence_path: Optional[str] = None
    qualifiers: Optional[Dict[str, List[Union[str, int, bool, float]]]] = None
    start: Optional[int] = None
    end: Optional[int] = None
    completely_within: Optional[bool] = None

    def to_annotation_collection(self, parent_or_seq_chunk_parent: Optional[Parent] = None) -> "AnnotationCollection":
        """Produce an :class:`~biocantor.gene.collections.AnnotationCollection` from this model."""
        if self.genes:
            genes = [gene.to_gene_interval(parent_or_seq_chunk_parent) for gene in self.genes]
        else:
            genes = None

        if self.feature_collections:
            feature_collections = [
                feat.to_feature_collection(parent_or_seq_chunk_parent) for feat in self.feature_collections
            ]
        else:
            feature_collections = None

        return AnnotationCollection(
            feature_collections=feature_collections,
            genes=genes,
            name=self.name,
            id=self.id,
            qualifiers=self.qualifiers,
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            start=self.start,
            end=self.end,
            completely_within=self.completely_within,
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    @staticmethod
    def from_annotation_collection(annotation_collection: AnnotationCollection) -> "AnnotationCollectionModel":
        """Convert back to :class:`~AnnotationCollectionModel`."""

        return AnnotationCollectionModel.Schema().load(annotation_collection.to_dict())
