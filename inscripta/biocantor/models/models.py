"""
Data models. These models allow for validation of inputs to a BioCantor model.
"""
from typing import List, Optional, ClassVar, Type
from uuid import UUID

from marshmallow import Schema
from marshmallow_dataclass import dataclass

from inscripta.biocantor.exc import InvalidCDSIntervalError, LocationException, ValidationException
from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds import CDSFrame, CDSInterval
from inscripta.biocantor.gene.collections import GeneInterval, FeatureIntervalCollection, AnnotationCollection
from inscripta.biocantor.gene.feature import FeatureInterval
from inscripta.biocantor.gene.transcript import TranscriptInterval
from inscripta.biocantor.location.location_impl import (
    CompoundInterval,
)
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import Parent


@dataclass
class BaseModel:
    """Base for all of the models."""

    Schema: ClassVar[Type[Schema]] = Schema  # noqa: F811

    class Meta:
        ordered = True


@dataclass
class FeatureIntervalModel(BaseModel):
    """Data model that allows construction of a :class:`~biocantor.gene.feature.FeatureInterval` object."""

    interval_starts: List[int]
    interval_ends: List[int]
    strand: Strand
    qualifiers: Optional[dict] = None
    sequence_symbol: Optional[str] = None
    sequence_guid: Optional[UUID] = None
    guid: Optional[UUID] = None
    feature_type: Optional[str] = None
    feature_symbol: Optional[str] = None
    feature_id: Optional[str] = None
    is_primary_feature: Optional[bool] = None

    def to_feature_interval(
        self,
        parent: Optional[Parent] = None,  # should have a sequence associated
    ) -> "FeatureInterval":
        """Construct from a :class:`~biocantor.models.FeatureIntervalModel`.

        A :class:`Parent` can be provided to allow the sequence-retrieval functions to work.
        """

        if len(self.interval_starts) != len(self.interval_ends):
            raise ValidationException("Number of interval starts does not match number of interval ends")

        location = CompoundInterval(
            self.interval_starts,
            self.interval_ends,
            self.strand,
            parent=parent,
        )

        return FeatureInterval(
            location=location,
            qualifiers=self.qualifiers,
            sequence_guid=self.sequence_guid,
            sequence_symbol=self.sequence_symbol,
            feature_type=self.feature_type,
            feature_symbol=self.feature_symbol,
            feature_id=self.feature_id,
            guid=self.guid,
            is_primary_feature=self.is_primary_feature,
        )

    @staticmethod
    def from_feature_interval(feature_interval: FeatureInterval) -> "FeatureIntervalModel":
        """Convert to a :class:`~biocantor.models.TranscriptIntervalModel`"""

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
    qualifiers: Optional[dict] = None
    is_primary_tx: Optional[bool] = None
    transcript_id: Optional[str] = None
    protein_id: Optional[str] = None
    transcript_symbol: Optional[str] = None
    transcript_type: Optional[Biotype] = None
    sequence_symbol: Optional[str] = None
    sequence_guid: Optional[UUID] = None
    guid: Optional[UUID] = None
    transcript_guid: Optional[UUID] = None

    def to_transcript_interval(
        self,
        parent: Optional[Parent] = None,  # should have a sequence associated
    ) -> "TranscriptInterval":
        """Construct from a :class:`~biocantor.models.TranscriptIntervalModel`.

        A :class:`Parent` can be provided to allow the sequence-retrieval functions to work.
        """

        if len(self.exon_starts) != len(self.exon_ends):
            raise LocationException("Number of exon starts does not match number of exon ends")

        location = CompoundInterval(self.exon_starts, self.exon_ends, self.strand, parent=parent)

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

            cds_interval = CompoundInterval(self.cds_starts, self.cds_ends, self.strand, parent=parent)
            if not cds_interval:
                raise InvalidCDSIntervalError("CDS must have a non-zero length")

            # length validation happens in the CDS constructor
            cds = CDSInterval(cds_interval, self.cds_frames)
        else:
            cds = None

        return TranscriptInterval(
            location=location,
            cds=cds,
            transcript_guid=self.guid,
            qualifiers=self.qualifiers,
            is_primary_tx=self.is_primary_tx,
            transcript_id=self.transcript_id,
            transcript_symbol=self.transcript_symbol,
            transcript_type=self.transcript_type if self.transcript_type else None,
            sequence_symbol=self.sequence_symbol,
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
    qualifiers: Optional[dict] = None
    sequence_symbol: Optional[str] = None
    sequence_guid: Optional[UUID] = None
    gene_guid: Optional[UUID] = None

    def to_gene_interval(self, parent: Optional[Parent] = None) -> GeneInterval:
        """Produce a :class:`GeneInterval` from a :class:`~biocantor.models.GeneIntervalModel`.

        This is the primary method of constructing a :class:`GeneInterval`.
        """

        transcripts = [tx.to_transcript_interval(parent) for tx in self.transcripts]

        return GeneInterval(
            transcripts=transcripts,
            guid=self.gene_guid,
            gene_id=self.gene_id,
            gene_symbol=self.gene_symbol,
            gene_type=self.gene_type,
            locus_tag=self.locus_tag,
            qualifiers=self.qualifiers,
            sequence_symbol=self.sequence_symbol,
            sequence_guid=self.sequence_guid,
            parent=parent,
        )

    @staticmethod
    def from_gene_interval(gene: GeneInterval) -> "GeneIntervalModel":
        return GeneIntervalModel.Schema().load(gene.to_dict())


@dataclass
class FeatureIntervalCollectionModel(BaseModel):
    """
    Data model that allows construction of :class:`~biocantor.gene.collections.FeatureCollection` object.

    This is a container for one or more :class:`~biocantor.gene.feature.FeatureInterval` objects.

    Feature types are arbitrary here, and not enumerated, because I am not trying
    to capture all of Sequence Ontology at this point.
    """

    feature_intervals: List[FeatureIntervalModel]
    feature_symbol: Optional[str] = None
    feature_id: Optional[str] = None
    feature_type: Optional[str] = None
    sequence_symbol: Optional[str] = None
    sequence_guid: Optional[UUID] = None
    feature_collection_guid: Optional[UUID] = None
    qualifiers: Optional[dict] = None

    def to_feature_collection(self, parent: Optional[Parent] = None) -> FeatureIntervalCollection:
        """Produce a feature collection from a :class:`~biocantor.models.FeatureIntervalCollectionModel`."""

        feature_intervals = [feat.to_feature_interval(parent) for feat in self.feature_intervals]

        return FeatureIntervalCollection(
            feature_intervals=feature_intervals,
            feature_symbol=self.feature_symbol,
            feature_id=self.feature_id,
            feature_type=self.feature_type,
            sequence_symbol=self.sequence_symbol,
            sequence_guid=self.sequence_guid,
            qualifiers=self.qualifiers,
            guid=self.feature_collection_guid,
            parent=parent,
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
    sequence_symbol: Optional[str] = None
    sequence_guid: Optional[UUID] = None
    sequence_path: Optional[str] = None
    qualifiers: Optional[dict] = None
    start: Optional[int] = None
    end: Optional[int] = None
    completely_within: Optional[bool] = None

    def to_annotation_collection(self, parent: Optional[Parent] = None) -> "AnnotationCollection":
        """Produce an :class:`AnnotationCollection` directly from lists of data."""
        if self.genes:
            genes = [gene.to_gene_interval(parent) for gene in self.genes]
        else:
            genes = None

        if self.feature_collections:
            feature_collections = [feat.to_feature_collection(parent) for feat in self.feature_collections]
        else:
            feature_collections = None

        return AnnotationCollection(
            feature_collections=feature_collections,
            genes=genes,
            name=self.name,
            qualifiers=self.qualifiers,
            sequence_symbol=self.sequence_symbol,
            sequence_guid=self.sequence_guid,
            start=self.start,
            end=self.end,
            completely_within=self.completely_within,
            parent=parent,
        )

    @staticmethod
    def from_annotation_collection(annotation_collection: AnnotationCollection) -> "AnnotationCollectionModel":
        """Convert back to :class:`AnnotationCollectionModel`."""

        return AnnotationCollectionModel.Schema().load(annotation_collection.to_dict())
