"""
This module contains :class:`~biocantor.gene.variants.VariantInterval`, which models diploid sequence variation.

This model is intended to be as simple as possible, and represent a single alternative haplotype. Variants
are always represented on the positive strand, and is loosely modeled after VCF files.
"""
from typing import Optional, Dict, Hashable, Any, Iterable, Set, List
from uuid import UUID

from inscripta.biocantor.exc import DuplicateFeatureError
from inscripta.biocantor.gene.interval import AbstractFeatureIntervalCollection, IntervalType
from inscripta.biocantor.gene.interval import AbstractInterval, AbstractFeatureInterval, QualifierValue
from inscripta.biocantor.io.bed import RGB, BED12
from inscripta.biocantor.io.gff3.rows import GFFRow
from inscripta.biocantor.location import Parent, SingleInterval, Strand
from inscripta.biocantor.sequence.sequence import Sequence, Alphabet
from inscripta.biocantor.util.bins import bins
from inscripta.biocantor.util.hashing import digest_object


class VariantInterval(AbstractFeatureInterval):
    _identifiers = ["variant_name", "variant_id"]

    def __init__(
        self,
        start: int,
        end: int,
        sequence: str,
        variant_type: str,
        phase_block: Optional[int],
        guid: Optional[UUID],
        variant_guid: Optional[UUID],
        variant_name: Optional[str],
        variant_id: Optional[str],
        qualifiers: Optional[Dict[Hashable, QualifierValue]] = None,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ):
        self._location = SingleInterval(start, end, Strand.PLUS, parent_or_seq_chunk_parent)
        self.sequence = Sequence(sequence, Alphabet.NT_STRICT_UNKNOWN)
        self.variant_type = variant_type
        self.phase_block = phase_block
        self.variant_name = variant_name
        self.variant_id = variant_id
        self.variant_guid = variant_guid
        # qualifiers come in as a List, convert to Set
        self._import_qualifiers_from_list(qualifiers)
        self.bin = bins(self.start, self.end, fmt="bed")

        if guid is None:
            self.guid = digest_object(
                start,
                end,
                self.qualifiers,
                self.sequence,
                self.variant_type,
                self.phase_block,
                self.variant_name,
                self.variant_guid,
            )
        else:
            self.guid = guid

    def __str__(self):
        return f"VariantInterval(({self.chromosome_location}), name={self.variant_name}, alt={self.sequence})"

    def __repr__(self):
        return "<{}>".format(str(self))

    def export_qualifiers(
        self, parent_qualifiers: Optional[Dict[Hashable, Set[Hashable]]] = None
    ) -> Dict[Hashable, Set[str]]:
        qualifiers = self._merge_qualifiers(parent_qualifiers)
        return qualifiers

    def to_bed12(
        self,
        score: Optional[int] = 0,
        rgb: Optional[RGB] = RGB(0, 0, 0),
        name: Optional[str] = "feature_name",
        chromosome_relative_coordinates: bool = True,
    ) -> BED12:
        raise NotImplementedError

    def to_gff(
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict] = None,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterable[GFFRow]:
        raise NotImplementedError

    def to_vcf(self):
        raise NotImplementedError

    def to_dict(self, chromosome_relative_coordinates: bool = True) -> Dict[str, Any]:
        if chromosome_relative_coordinates:
            start = self.start
            end = self.end
        else:
            b = list(self.relative_blocks)[0]
            start = b.start
            end = b.end
        return dict(
            start=start,
            end=end,
            sequence=str(self.sequence),
            variant_type=self.variant_type,
            phase_block=self.phase_block,
            guid=self.guid,
            variant_guid=self.variant_guid,
            variant_name=self.variant_name,
            variant_id=self.variant_id,
            qualifiers=self._export_qualifiers_to_list(),
        )

    @staticmethod
    def from_dict(vals: Dict[str, Any], parent_or_seq_chunk_parent: Optional[Parent] = None) -> "VariantInterval":
        return VariantInterval(
            vals["start"],
            vals["end"],
            Sequence(vals["sequence"], Alphabet.NT_STRICT_UNKNOWN),
            vals["variant_type"],
            vals["phase_block"],
            vals["guid"],
            vals["variant_guid"],
            vals["variant_name"],
            vals["variant_id"],
            vals["qualifiers"],
        )

    @property
    def id(self) -> str:
        """Returns the ID of this feature. Provides a shared API across genes/transcripts and features."""
        return self.variant_id

    @property
    def name(self) -> str:
        """Returns the name of this feature. Provides a shared API across genes/transcripts and features."""
        return self.variant_name


class VariantIntervalCollection(AbstractFeatureIntervalCollection):
    """
    A container for many :class:`VariantInterval`. Assumes that the variants are all on the same haplotype.
    """

    interval_type = IntervalType.VARIANT
    _identifiers = ["variant_collection_name", "variant_collection_id"]

    def __init__(
        self,
        variant_intervals: List[VariantInterval],
        variant_collection_name: Optional[str] = None,
        variant_collection_id: Optional[str] = None,
        sequence_name: Optional[str] = None,
        sequence_guid: Optional[UUID] = None,
        guid: Optional[UUID] = None,
        qualifiers: Optional[Dict[Hashable, List[QualifierValue]]] = None,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ):
        self.variant_intervals = variant_intervals
        self.variant_collection_name = variant_collection_name
        self.variant_collection_id = variant_collection_id
        self.sequence_name = sequence_name
        self.sequence_guid = sequence_guid
        # qualifiers come in as a List, convert to Set
        self._import_qualifiers_from_list(qualifiers)
        self.start = self.genomic_start = min(f.start for f in self.variant_intervals)
        self.end = self.genomic_end = max(f.end for f in self.variant_intervals)
        self._initialize_location(self.start, self.end, parent_or_seq_chunk_parent)

        self.variant_types = {x.variant_type for x in self.variant_intervals}

        if guid is None:
            self.guid = digest_object(
                self.chunk_relative_location,
                self.variant_collection_name,
                self.variant_collection_id,
                self.sequence_name,
                self.qualifiers,
                self.children_guids,
            )
        else:
            self.guid = guid

        self.guid_map = {}
        for feat in self.variant_intervals:
            if feat.guid in self.guid_map:
                raise DuplicateFeatureError(f"Guid {feat.guid} found more than once in this VariantIntervalCollection")
            self.guid_map[feat.guid] = feat

    def iter_children(self) -> Iterable["AbstractInterval"]:
        yield from self.variant_intervals

    def children_guids(self) -> Set[UUID]:
        return set(self.guid_map.keys())

    def query_by_guids(self, ids: List[UUID]) -> "VariantIntervalCollection":
        variant_intervals = [self.guid_map[i] for i in ids if i in self.guid_map]
        if variant_intervals:
            return VariantIntervalCollection(
                variant_intervals=variant_intervals,
                variant_collection_name=self.variant_collection_name,
                variant_collection_id=self.variant_collection_id,
                qualifiers=self._export_qualifiers_to_list(),
                sequence_name=self.sequence_name,
                sequence_guid=self.sequence_guid,
                guid=self.guid,
                parent_or_seq_chunk_parent=self.chunk_relative_location.parent,
            )

    def to_dict(self, chromosome_relative_coordinates: bool = True) -> Dict[str, Any]:
        """Convert to a dict usable by :class:`~biocantor.io.models.VariantIntervalCollectionModel`."""
        return dict(
            variant_intervals=[var.to_dict(chromosome_relative_coordinates) for var in self.variant_intervals],
            variant_collection_name=self.variant_collection_name,
            variant_collection_id=self.variant_collection_id,
            qualifiers=self._export_qualifiers_to_list(),
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            variant_collection_guid=self.guid,
        )

    @staticmethod
    def from_dict(
        vals: Dict[str, Any], parent_or_seq_chunk_parent: Optional[Parent] = None
    ) -> "VariantIntervalCollection":
        """Build an :class:`VariantIntervalCollection` from a dictionary representation"""
        return VariantIntervalCollection(
            variant_intervals=[
                VariantInterval.from_dict(x, parent_or_seq_chunk_parent) for x in vals["variant_intervals"]
            ],
            variant_collection_name=vals["variant_collection_name"],
            variant_collection_id=vals["variant_collection_id"],
            qualifiers=vals["qualifiers"],
            sequence_name=vals["sequence_name"],
            sequence_guid=vals["sequence_guid"],
            guid=vals["feature_collection_guid"],
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    def to_gff(self, chromosome_relative_coordinates: bool = True) -> Iterable[GFFRow]:
        raise NotImplementedError

    @property
    def id(self) -> str:
        return self.variant_collection_id

    @property
    def name(self) -> str:
        return self.variant_collection_name
