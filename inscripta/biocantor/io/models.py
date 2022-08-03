"""
Data models. These models allow for validation of inputs to a BioCantor model, acting as a JSON schema for serializing
and deserializing the models.
"""
from typing import List, Optional, ClassVar, Type, Dict, Union
from dataclasses import field
from uuid import UUID

from marshmallow import Schema, post_dump  # noqa: F401
from marshmallow_dataclass import dataclass

from inscripta.biocantor.gene import GeneInterval, FeatureIntervalCollection
from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds_frame import CDSFrame
from inscripta.biocantor.gene.collections import AnnotationCollection
from inscripta.biocantor.gene.feature import FeatureInterval
from inscripta.biocantor.gene.transcript import TranscriptInterval
from inscripta.biocantor.gene.variants import VariantInterval, VariantIntervalCollection
from inscripta.biocantor.io.exc import InvalidInputError
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import Parent
from inscripta.biocantor.sequence.sequence import Alphabet, SequenceType


@dataclass
class BaseModel:
    """Base for all of the models."""

    Schema: ClassVar[Type[Schema]] = Schema  # noqa: F811

    class Meta:
        ordered = True


@dataclass
class ParentModel(BaseModel):
    """Data model that allows construction of a :class:`~biocantor.gene.parent.Parent` object."""

    seq: Optional[str] = None
    alphabet: Optional[Alphabet] = Alphabet.NT_EXTENDED_GAPPED
    sequence_name: Optional[str] = None
    type: Optional[Union[SequenceType, str]] = SequenceType.CHROMOSOME
    start: Optional[int] = None
    end: Optional[int] = None
    strand: Optional[Strand] = None

    def to_parent(self) -> Parent:
        if self.seq is None:
            return Parent(sequence_type=self.type, id=self.sequence_name)

        # avoid circular imports
        from inscripta.biocantor.io.parser import seq_chunk_to_parent, seq_to_parent

        seq_type = SequenceType.sequence_type_str_to_type(self.type)
        if seq_type == SequenceType.SEQUENCE_CHUNK:
            if not self.sequence_name:
                raise InvalidInputError("Cannot construct sequence chunk parent without a sequence name")
            elif self.start is None or self.end is None:
                raise InvalidInputError("Cannot construct sequence chunk parent without chunk start/end positions")
            # the need for fn_call and fn here is to allow seq_to_parent or seq_chunk_to_parent to retain
            # their default parameters when None is value here
            fn = seq_chunk_to_parent
            fn_call = dict(
                seq=self.seq,
                sequence_name=self.sequence_name,
                start=self.start,
                end=self.end,
                strand=self.strand,
                alphabet=self.alphabet,
            )
        else:
            fn = seq_to_parent
            fn_call = dict(seq=self.seq, alphabet=self.alphabet, seq_id=self.sequence_name, seq_type=seq_type)
        return fn(**{k: v for k, v in fn_call.items() if v is not None})


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
        return FeatureInterval(
            self.interval_starts,
            self.interval_ends,
            self.strand,
            qualifiers=self.qualifiers,
            sequence_guid=self.sequence_guid,
            sequence_name=self.sequence_name,
            feature_types=self.feature_types,
            feature_name=self.feature_name,
            feature_id=self.feature_id,
            guid=self.feature_interval_guid,
            feature_guid=self.feature_guid,
            is_primary_feature=self.is_primary_feature,
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
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
    product: Optional[str] = None
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

        return TranscriptInterval(
            exon_starts=self.exon_starts,
            exon_ends=self.exon_ends,
            strand=self.strand,
            cds_starts=self.cds_starts,
            cds_ends=self.cds_ends,
            cds_frames=self.cds_frames,
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
            product=self.product,
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    @staticmethod
    def from_transcript_interval(transcript_interval: TranscriptInterval) -> "TranscriptIntervalModel":
        """Convert to a :class:`~biocantor.io.models.TranscriptIntervalModel`"""

        return TranscriptIntervalModel.Schema().load(transcript_interval.to_dict())


@dataclass
class VariantIntervalModel(BaseModel):
    """
    Data model that allows construction of :class:`~biocantor.gene.variants.VariantInterval` object.
    """

    start: int
    end: int
    sequence: str
    variant_type: str
    phase_block: Optional[int] = None
    variant_interval_guid: Optional[UUID] = None
    variant_guid: Optional[UUID] = None
    variant_name: Optional[str] = None
    variant_id: Optional[str] = None
    qualifiers: Optional[Dict[str, List[Union[str, int, bool, float]]]] = None

    @staticmethod
    def from_variant_interval(variant_interval: VariantInterval) -> "VariantIntervalModel":
        """Convert to a :class:`~biocantor.io.models.VariantIntervalModel`"""

        return VariantIntervalModel.Schema().load(variant_interval.to_dict())

    def to_variant_interval(self, parent_or_seq_chunk_parent: Optional[Parent] = None) -> VariantInterval:
        return VariantInterval(
            self.start,
            self.end,
            self.sequence,
            self.variant_type,
            self.phase_block,
            self.variant_interval_guid,
            self.variant_guid,
            self.variant_name,
            self.variant_id,
            self.qualifiers,
            parent_or_seq_chunk_parent,
        )


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
class VariantIntervalCollectionModel(BaseModel):
    """
    Data model that allows construction of :class:`~biocantor.gene.collections.VariantCollection` object.

    This is a container for one or more :class:`~biocantor.gene.variants.VariantInterval` objects.

    VariantIntervalCollection stores one or more variants on a phase block.
    """

    variant_intervals: List[VariantIntervalModel]
    variant_collection_name: Optional[str] = None
    variant_collection_id: Optional[str] = None
    sequence_name: Optional[str] = None
    sequence_guid: Optional[UUID] = None
    variant_collection_guid: Optional[UUID] = None
    qualifiers: Optional[Dict[str, List[Union[str, int, bool, float]]]] = None

    @staticmethod
    def from_variant_interval_collection(
        variant_collection: VariantIntervalCollection,
    ) -> "VariantIntervalCollectionModel":
        """Convert to a :class:`~biocantor.io.models.VariantIntervalModel`"""

        return VariantIntervalCollectionModel.Schema().load(variant_collection.to_dict())

    def to_variant_interval_collection(
        self, parent_or_seq_chunk_parent: Optional[Parent] = None
    ) -> VariantIntervalCollection:
        """Produce a variant collection from a :class:`VariantIntervalCollectionModel`."""

        variant_intervals = [
            variant.to_variant_interval(parent_or_seq_chunk_parent) for variant in self.variant_intervals
        ]
        return VariantIntervalCollection(
            variant_intervals,
            self.variant_collection_name,
            self.variant_collection_id,
            self.sequence_name,
            self.sequence_guid,
            self.variant_collection_guid,
            self.qualifiers,
            parent_or_seq_chunk_parent,
        )


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

    feature_collections: List[FeatureIntervalCollectionModel] = field(default_factory=list)
    genes: List[GeneIntervalModel] = field(default_factory=list)
    variant_collections: List[VariantIntervalCollectionModel] = field(default_factory=list)
    name: Optional[str] = None
    id: Optional[str] = None
    sequence_name: Optional[str] = None
    sequence_guid: Optional[UUID] = None
    sequence_path: Optional[str] = None
    qualifiers: Optional[Dict[str, List[Union[str, int, bool, float]]]] = None
    start: Optional[int] = None
    end: Optional[int] = None
    completely_within: Optional[bool] = None
    parent_or_seq_chunk_parent: Optional[ParentModel] = None

    @post_dump(pass_original=True, pass_many=False)
    def post_dump(self, data, model, many=False):
        """
        If the object being dumped is an AnnotationCollection, convert the ``_parent_or_seq_chunk_parent``
        value to ``parent_or_seq_chunk_parent``.
        """
        if isinstance(model, AnnotationCollection):
            data["parent_or_seq_chunk_parent"] = model._parent_to_dict()
        return data

    def to_annotation_collection(self, parent_or_seq_chunk_parent: Optional[Parent] = None) -> "AnnotationCollection":
        """Produce an :class:`~biocantor.gene.collections.AnnotationCollection` from this model."""

        if not parent_or_seq_chunk_parent and self.parent_or_seq_chunk_parent:
            parent_or_seq_chunk_parent = self.parent_or_seq_chunk_parent.to_parent()

        genes = [gene.to_gene_interval(parent_or_seq_chunk_parent) for gene in self.genes]
        feature_collections = [
            feat.to_feature_collection(parent_or_seq_chunk_parent) for feat in self.feature_collections
        ]
        variant_collections = [
            var.to_variant_interval_collection(parent_or_seq_chunk_parent) for var in self.variant_collections
        ]

        return AnnotationCollection(
            feature_collections=feature_collections,
            genes=genes,
            variant_collections=variant_collections,
            name=self.name,
            id=self.id,
            qualifiers=self.qualifiers,
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            sequence_path=self.sequence_path,
            start=self.start,
            end=self.end,
            completely_within=self.completely_within,
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    @staticmethod
    def from_annotation_collection(
        annotation_collection: AnnotationCollection,
        chromosome_relative_coordinates: bool = True,
        export_parent: bool = False,
    ) -> "AnnotationCollectionModel":
        """Convert back to :class:`~AnnotationCollectionModel`."""

        return AnnotationCollectionModel.Schema().load(
            annotation_collection.to_dict(
                chromosome_relative_coordinates=chromosome_relative_coordinates, export_parent=export_parent
            )
        )
