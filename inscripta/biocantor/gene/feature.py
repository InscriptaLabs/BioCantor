"""
Object representation of features. Includes an abstract feature class that is also used by transcripts.

Each object is capable of exporting itself to BED and GFF3.
"""
from abc import ABC, abstractmethod
from functools import lru_cache
from typing import Optional, Any, Union, Dict, List, Set, Iterable, Hashable
from uuid import UUID

from inscripta.biocantor.exc import (
    EmptyLocationException,
)
from inscripta.biocantor.gene.cds import CDSPhase
from inscripta.biocantor.io.bed import BED12, RGB
from inscripta.biocantor.io.gff3.constants import GFF_SOURCE, NULL_COLUMN, BioCantorFeatureTypes, BioCantorQualifiers
from inscripta.biocantor.io.gff3.rows import GFFAttributes, GFFRow
from inscripta.biocantor.location.location import Location
from inscripta.biocantor.location.location_impl import SingleInterval
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent.parent import Parent
from inscripta.biocantor.sequence.sequence import Sequence
from inscripta.biocantor.util.bins import bins
from inscripta.biocantor.util.hashing import digest_object
from inscripta.biocantor.util.object_validation import ObjectValidation


class AbstractFeatureInterval(ABC):
    """This is a wrapper over :class:`~biocantor.location.Location` that adds metadata coordinate transformation
    QOL functions."""

    location: Location
    _identifiers: List[Union[str, UUID]]
    qualifiers: Dict[Hashable, Set[Hashable]]
    guid: Optional[UUID] = None
    sequence_guid: Optional[UUID] = None
    sequence_name: Optional[str] = None
    _is_primary_feature: Optional[bool] = None

    @abstractmethod
    def export_qualifiers(
        self, parent_qualifiers: Optional[Dict[Hashable, Set[Hashable]]] = None
    ) -> Dict[Hashable, Set[Hashable]]:
        """Exports qualifiers for GFF3 or GenBank export. This merges top level keys with the arbitrary values"""

    @property
    def start(self) -> int:
        """Returns start position."""
        return self.location.start

    @property
    def end(self) -> int:
        """Returns end position."""
        return self.location.end

    @property
    def strand(self) -> Strand:
        """Returns strand of location."""
        return self.location.strand

    @property
    def identifiers(self) -> Set[Union[str, UUID]]:
        """Returns the identifiers for this FeatureInterval, if they exist"""
        return {getattr(self, i) for i in self._identifiers if getattr(self, i) is not None}

    @property
    def identifiers_dict(self) -> Dict[str, Union[str, UUID]]:
        """Returns the identifiers and their keys for this FeatureInterval, if they exist"""
        return {key: getattr(self, key) for key in self._identifiers if getattr(self, key) is not None}

    @property
    def is_primary_feature(self) -> bool:
        """Is this the primary feature?"""
        return self._is_primary_feature is True

    def __len__(self):
        return len(self.location)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.__dict__ == other.__dict__

    def __hash__(self):
        """Produces a hash, which is the GUID."""
        return hash(self.guid)

    @abstractmethod
    def update_parent(self, parent: Parent):
        """Update parent"""

    @abstractmethod
    def to_bed12(
        self, score: Optional[int] = 0, rgb: Optional[RGB] = RGB(0, 0, 0), name: Optional[str] = "feature_name"
    ) -> BED12:
        """Write a BED12 format representation of this :class:`AbstractFeatureInterval`.

        Both of these optional arguments are specific to the BED12 format.

        Args:
            score: An optional score associated with a interval. UCSC requires an integer between 0 and 1000.
            rgb: An optional RGB string for visualization on a browser. This allows you to have multiple colors
                on a single UCSC track.
            name: Which identifier in this record to use as 'name'. feature_name to guid. If the supplied string
                is not a valid attribute, it is used directly.

        Return:
            A :class:`~biocantor.io.bed.BED12` object.
        """

    @abstractmethod
    def to_gff(self, parent: Optional[str] = None, parent_qualifiers: Optional[Dict] = None) -> Iterable[GFFRow]:
        """Writes a GFF format list of lists for this feature.

        The additional qualifiers are used when writing a hierarchical relationship back to files. GFF files
        are easier to work with if the children features have the qualifiers of their parents.

        Args:
            parent: ID of the Parent of this transcript.
            parent_qualifiers: Directly pull qualifiers in from this dictionary.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`
        """

    def sequence_pos_to_feature(self, pos: int) -> int:
        """Converts sequence position to relative position along this feature."""
        return self.location.parent_to_relative_pos(pos)

    def sequence_interval_to_feature(self, chr_start: int, chr_end: int, chr_strand: Strand) -> Location:
        """Converts a contiguous interval on the sequence to a relative location within this feature."""
        return self.location.parent_to_relative_location(SingleInterval(chr_start, chr_end, chr_strand))

    def feature_pos_to_sequence(self, pos: int) -> int:
        """Converts a relative position along this feature to sequence coordinate."""
        return self.location.relative_to_parent_pos(pos)

    def feature_interval_to_sequence(self, rel_start: int, rel_end: int, rel_strand: Strand) -> Location:
        """Converts a contiguous interval relative to this feature to a spliced location on the sequence."""
        return self.location.relative_interval_to_parent_location(rel_start, rel_end, rel_strand)

    @lru_cache(maxsize=1)
    def get_spliced_sequence(self) -> Sequence:
        """Returns the feature's *spliced*, *stranded* sequence."""
        return self.location.extract_sequence()

    @lru_cache(maxsize=1)
    def get_reference_sequence(self) -> Sequence:
        """Returns the feature's *unspliced*, *positive strand* genomic sequence."""
        ObjectValidation.require_location_has_parent_with_sequence(self.location)
        return self.location.parent.sequence[self.location.start : self.location.end]

    @lru_cache(maxsize=1)
    def get_genomic_sequence(self) -> Sequence:
        """Returns the feature's *unspliced*, *stranded* (transcription orientation) genomic sequence."""
        seq = self.location.parent.sequence[self.location.start : self.location.end]
        if self.strand == Strand.PLUS:
            return seq
        else:
            return seq.reverse_complement()

    def _merge_qualifiers(
        self, other_qualifiers: Optional[Dict[Hashable, Set[Hashable]]] = None
    ) -> Dict[Hashable, Set[Hashable]]:
        """Merges this Interval's qualifiers dictionary with a new one, removing redundancy."""
        merged = self.qualifiers.copy()
        if other_qualifiers:
            for key, vals in other_qualifiers.items():
                if key not in merged:
                    merged[key] = set()
                merged[key].update(vals)
        return merged


class FeatureInterval(AbstractFeatureInterval):
    """FeatureIntervals are generic intervals. These can be used to model genome promoters,
    open chromatin sites, etc.
    """

    _identifiers = ["feature_name", "feature_id"]

    def __init__(
        self,
        location: Location,  # exons
        qualifiers: Optional[Dict[Hashable, Set[Hashable]]] = None,
        sequence_guid: Optional[UUID] = None,
        sequence_name: Optional[str] = None,
        feature_type: Optional[str] = None,
        feature_name: Optional[str] = None,
        feature_id: Optional[str] = None,
        guid: Optional[UUID] = None,
        is_primary_feature: Optional[bool] = None,
    ):
        self.location = location  # genomic CompoundInterval
        self.sequence_guid = sequence_guid
        self.sequence_name = sequence_name
        self.feature_type = feature_type
        self.feature_name = feature_name
        self.feature_id = feature_id

        if qualifiers:
            self.qualifiers = qualifiers
        else:
            self.qualifiers = {}

        self.bin = bins(self.start, self.end, fmt="bed")
        self._is_primary_feature = is_primary_feature

        if guid is None:
            self.guid = digest_object(
                self.location,
                self.qualifiers,
                self.sequence_name,
                self.feature_type,
                self.feature_name,
                self.feature_id,
                self.is_primary_feature,
            )
        else:
            self.guid = guid

        if self.location.parent:
            ObjectValidation.require_location_has_parent_with_sequence(self.location)

    def __str__(self):
        return f"FeatureInterval(({self.location}), name={self.feature_name})"

    def __repr__(self):
        return "<{}>".format(str(self))

    def to_dict(self) -> Dict[str, Any]:
        """Convert to a dict usable by :class:`FeatureIntervalModel`."""
        exon_starts, exon_ends = list(zip(*([x.start, x.end] for x in self.location.blocks)))
        return dict(
            interval_starts=exon_starts,
            interval_ends=exon_ends,
            strand=self.strand.name,
            qualifiers=self.qualifiers if self.qualifiers else None,
            feature_id=self.feature_id,
            feature_name=self.feature_name,
            feature_type=self.feature_type,
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            feature_interval_guid=self.guid,
            is_primary_feature=self._is_primary_feature,
        )

    def update_parent(self, parent: Parent):
        """Change parent to this new parent."""
        self.location = self.location.reset_parent(parent)

    def intersect(
        self,
        location: Location,
        new_guid: Optional[UUID] = None,
        new_qualifiers: Optional[dict] = None,
    ) -> "FeatureInterval":
        """Returns a new FeatureInterval representing the intersection of this FeatureInterval's location with the
        other location.

        Strand of the other location is ignored; returned FeatureInterval is on the same strand as this FeatureInterval.
        """
        if not new_qualifiers:
            new_qualifiers = self.qualifiers

        location_same_strand = location.reset_strand(self.location.strand)
        intersection = self.location.intersection(location_same_strand)

        if intersection.is_empty:
            raise EmptyLocationException("Can't intersect disjoint intervals")

        return FeatureInterval(location=intersection, guid=new_guid, qualifiers=new_qualifiers)

    def export_qualifiers(
        self, parent_qualifiers: Optional[Dict[Hashable, Set[Hashable]]] = None
    ) -> Dict[Hashable, Set[Hashable]]:
        """Exports qualifiers for GFF3/GenBank export"""
        qualifiers = self._merge_qualifiers(parent_qualifiers)
        for key, val in [
            [BioCantorQualifiers.FEATURE_SYMBOL.value, self.feature_name],
            [BioCantorQualifiers.FEATURE_ID.value, self.feature_id],
            [BioCantorQualifiers.FEATURE_TYPE.value, self.feature_type],
        ]:
            if not val:
                continue
            if key not in qualifiers:
                qualifiers[key] = set()
            qualifiers[key].add(val)
        return qualifiers

    def to_gff(
        self, parent: Optional[str] = None, parent_qualifiers: Optional[Dict[Hashable, Set[Hashable]]] = None
    ) -> Iterable[GFFRow]:
        """Writes a GFF format list of lists for this feature.

        The additional qualifiers are used when writing a hierarchical relationship back to files. GFF files
        are easier to work with if the children features have the qualifiers of their parents.

        Args:
            parent: ID of the Parent of this transcript.
            parent_qualifiers: Directly pull qualifiers in from this dictionary.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`
        """
        qualifiers = self.export_qualifiers(parent_qualifiers)

        feature_id = str(self.guid) if self.guid else str(digest_object(self))

        attributes = GFFAttributes(id=feature_id, qualifiers=qualifiers, name=self.feature_name, parent=parent)

        # transcript feature
        row = GFFRow(
            self.sequence_name,
            GFF_SOURCE,
            BioCantorFeatureTypes.TRANSCRIPT,
            self.start + 1,
            self.end,
            NULL_COLUMN,
            self.strand,
            CDSPhase.NONE,
            attributes,
        )
        yield row

        # start adding exon features
        # re-use qualifiers, updating ID each time
        for i, block in enumerate(self.location.blocks, 1):
            attributes = GFFAttributes(
                id=f"exon-{feature_id}-{i}", qualifiers=qualifiers, name=self.feature_name, parent=feature_id
            )
            row = GFFRow(
                self.sequence_name,
                GFF_SOURCE,
                BioCantorFeatureTypes.EXON,
                block.start + 1,
                block.end,
                NULL_COLUMN,
                self.strand,
                CDSPhase.NONE,
                attributes,
            )
            yield row

    def to_bed12(
        self, score: Optional[int] = 0, rgb: Optional[RGB] = RGB(0, 0, 0), name: Optional[str] = "feature_name"
    ) -> BED12:
        block_sizes = [b.end - b.start for b in self.location.blocks]
        block_starts = [b.start - self.start for b in self.location.blocks]
        return BED12(
            self.sequence_name,
            self.start,
            self.end,
            getattr(self, name, name),
            score,
            self.strand,
            0,  # thickStart always 0 for non-coding
            0,  # thickEnd always 0 for non-coding
            rgb,
            len(self.location.blocks),
            block_sizes,
            block_starts,
        )
