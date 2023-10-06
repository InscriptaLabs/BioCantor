"""
This module contains :class:`~biocantor.gene.variants.VariantInterval`, which models diploid sequence variation.

This model is intended to be as simple as possible, and represent a single alternative haplotype. Variants
are always represented on the positive strand, and is loosely modeled after VCF files.

Variants must be represented by intervals of at least length 1 in all cases. It is better to left-pad indels,
but it is not required.

Here are examples of how variants can be represented.

A SNV:

.. code-block:: python

    # ref: ACTCTCTCTATCTCATCCAC
    # alt: AGTCTCTCTATCTCATCCAC
    snp_1 = VariantInterval(start=1, end=2, sequence="G", variant_type="SNV")

A GG insertion at position 5

.. code-block:: python

    # ref: ACTCT  CTCTATCTCATCCAC
    # alt: ACTCTGGCTCTATCTCATCCAC
    insertion_5 = VariantInterval(start=5, end=6, sequence="GGC", variant_type="insertion")

A left-padded deletion from 10 to 13:

.. code-block:: python

    # ref: ACTCTCTCTATCTCATCCAC
    # alt: ACTCTCTCTAT  CATCCAC
    deletion_11_13 = VariantInterval(start=10, end=13, sequence="T", variant_type="deletion")

A deletion from 13 to 15 without padding:

.. code-block:: python

    # ref: ACTCTCTCTATCTCATCCAC
    # alt: ACTCTCTCTATCT  TCCAC
    deletion_13_15 = VariantInterval(start=13, end=15, sequence="", variant_type="deletion")

"""
from typing import Optional, Dict, Hashable, Any, Iterable, Iterator, Set, List, Union
from uuid import UUID

from biocantor.exc import (
    DuplicateFeatureError,
    LocationOverlapException,
    EmptyLocationException,
    NullSequenceException,
)
from biocantor.gene.interval import (
    AbstractFeatureIntervalCollection,
    IntervalType,
    AbstractInterval,
    AbstractFeatureInterval,
    QualifierValue,
)
from biocantor.io.bed import RGB, BED12
from biocantor.io.gff3.rows import GFFRow
from biocantor.location import Parent, SingleInterval, Strand, Location, CompoundInterval, EmptyLocation
from biocantor.sequence.sequence import Sequence, Alphabet, SequenceType
from biocantor.util.bins import bins
from biocantor.util.hashing import digest_object


class VariantInterval(AbstractFeatureInterval):
    _identifiers = ["variant_name", "variant_id"]

    def __init__(
        self,
        start: int,
        end: int,
        sequence: str,
        variant_type: str,
        phase_block: Optional[int] = None,
        guid: Optional[UUID] = None,
        variant_guid: Optional[UUID] = None,
        variant_name: Optional[str] = None,
        variant_id: Optional[str] = None,
        qualifiers: Optional[Dict[Hashable, QualifierValue]] = None,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ):
        if start == end:
            raise EmptyLocationException(
                "Variants must be defined over a window of at least 1bp. Indels should be left-padded."
            )
        self._location = VariantInterval.initialize_location([start], [end], Strand.PLUS, parent_or_seq_chunk_parent)
        self._parent_or_seq_chunk_parent = parent_or_seq_chunk_parent
        self.sequence = Sequence(sequence, Alphabet.NT_STRICT_UNKNOWN)
        self.variant_type = variant_type
        self.phase_block = phase_block
        self.variant_name = variant_name
        self.variant_id = variant_id
        self.variant_guid = variant_guid
        # qualifiers come in as a List, convert to Set
        self._import_qualifiers_from_list(qualifiers)
        self.bin = bins(start, end, fmt="bed")
        self.start = self.genomic_start = start
        self.end = self.genomic_end = end
        self._genomic_starts = [start]
        self._genomic_ends = [end]
        # variants are always + strand
        self._strand = Strand.PLUS

        # lazy loaded values
        self._alternative_sequence = None
        self._parent_with_alternative_sequence = None

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

    def to_gtf(
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict] = None,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterator[GFFRow]:
        raise NotImplementedError

    def to_gff(
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict] = None,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterator[GFFRow]:
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
            vals["sequence"],
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

    @property
    def alternative_genomic_sequence(self) -> Sequence:
        """Edited version of the original genomic sequence"""
        if not self.has_sequence:
            raise NullSequenceException("This VariantInterval has no sequence information")
        if not self._alternative_sequence:
            original_sequence = self.chunk_relative_location.parent.sequence
            original_sequence_str = str(original_sequence)
            new_seq_string = "".join(
                (
                    original_sequence_str[: self.chunk_relative_location.start],
                    str(self.sequence),
                    original_sequence_str[self.chunk_relative_location.end :],
                )
            )
            self._alternative_sequence = Sequence(
                data=new_seq_string,
                alphabet=original_sequence.alphabet,
                id=None,
                type=self.variant_type,
                parent=None,
                validate_alphabet=False,
            )
        return self._alternative_sequence

    @property
    def parent_with_alternative_sequence(self) -> Parent:
        if not self.has_sequence:
            raise NullSequenceException("This VariantInterval has no sequence information")
        if self._parent_with_alternative_sequence is None:
            # have to import here to avoid circular imports
            from biocantor.io.parser import seq_chunk_to_parent, seq_to_parent

            if self.chunk_relative_location.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
                self._parent_with_alternative_sequence = seq_chunk_to_parent(
                    str(self.alternative_genomic_sequence),
                    self.chunk_relative_location.parent_id,
                    self.chunk_relative_location.parent.parent.location.start,
                    self.chunk_relative_location.parent.parent.location.start + len(self.alternative_genomic_sequence),
                )
            else:
                self._parent_with_alternative_sequence = seq_to_parent(str(self.alternative_genomic_sequence))
        return self._parent_with_alternative_sequence

    @property
    def length_difference(self):
        return len(self.sequence) - len(self.chromosome_location)

    def lift_over_location(self, location: Location) -> Location:
        """
        Construct a new Location that takes the alternative sequence defined by this VariantInterval into account.

        The Location can be chunk-relative or chromosome-relative. It will be returned relative to the coordinate
        system of this VariantInterval.
        """
        if location is EmptyLocation():
            return EmptyLocation()

        # all liftovers happen in chromosome coordinates
        if location.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            location = location.lift_over_to_first_ancestor_of_type(SequenceType.CHROMOSOME)

        if len(self.chromosome_location) == len(self.sequence) or location.end <= self.chromosome_location.start:
            return self.liftover_location_to_seq_chunk_parent(location, self.parent_with_alternative_sequence)

        if type(location) is SingleInterval:
            new_loc = self._lift_over_chromosome_location_single_interval(location)
        elif type(location) is CompoundInterval:
            new_loc = self._lift_over_chromosome_location_compound_interval(location)
        else:
            raise NotImplementedError("Location type {} not supported".format(str(type(location))))
        # this lifts the chromosome coordinates back onto chunk coordinates, if we are chunk-relative
        if self.has_sequence:
            return self.liftover_location_to_seq_chunk_parent(new_loc, self.parent_with_alternative_sequence)
        else:
            return new_loc

    def _lift_over_chromosome_location_single_interval(self, location: SingleInterval) -> Location:
        len_diff = self.length_difference
        old_start = location.start
        old_end = location.end
        if len_diff >= 0:
            new_start = old_start if old_start < self.chromosome_location.end else old_start + len_diff
            new_end = old_end if old_end < self.chromosome_location.end else old_end + len_diff
        else:
            if (
                self.chromosome_location.start + len(self.sequence)
                <= old_start
                <= old_end
                <= self.chromosome_location.end
            ):
                return EmptyLocation()
            len_left_side_deleted = (
                self.chromosome_location.end - old_start
                if self.chromosome_location.start + len(self.sequence) <= old_start < self.chromosome_location.end
                else 0
            )
            new_start = (
                old_start
                if old_start < self.chromosome_location.end + len_diff
                else old_start + len_diff + len_left_side_deleted
            )
            new_end = (
                old_end
                if old_end <= self.chromosome_location.end + len_diff
                else old_end + max(len_diff, len_diff - old_end + self.chromosome_location.end)
            )
        return SingleInterval(new_start, new_end, location.strand)

    def _lift_over_chromosome_location_compound_interval(self, location: CompoundInterval) -> Location:
        lifted_single_intervals = []
        for single_interval in location.blocks:
            try:
                lifted_interval = self._lift_over_chromosome_location_single_interval(single_interval)
                if lifted_interval is not EmptyLocation():
                    lifted_single_intervals.append(lifted_interval)
            except EmptyLocationException:
                continue
        if lifted_single_intervals:
            return (
                lifted_single_intervals[0]
                if len(lifted_single_intervals) == 1
                else CompoundInterval.from_single_intervals(lifted_single_intervals).optimize_blocks()
            )
        else:
            return EmptyLocation()


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
        self.variant_intervals = sorted(variant_intervals, key=lambda x: x.start)
        # validate the variant intervals for not being overlapping
        for i in range(len(self.variant_intervals) - 1):
            if self.variant_intervals[i].chunk_relative_location.has_overlap(
                self.variant_intervals[i + 1].chunk_relative_location
            ):
                raise LocationOverlapException("VariantInterval within a VariantIntervalCollection must not overlap")

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
        self._parent_or_seq_chunk_parent = parent_or_seq_chunk_parent

        # lazy-loaded property
        self._alternative_genomic_sequence = None
        self._parent_with_alternative_sequence = None

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

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(identifiers={self.identifiers}, "
            f"Intervals:{','.join(str(f) for f in self.variant_intervals)})"
        )

    def iter_children(self) -> Iterable["AbstractInterval"]:
        yield from self.variant_intervals

    @property
    def children_guids(self) -> Set[UUID]:
        return {x.guid for x in self.variant_intervals}

    def query_by_guids(self, id_or_ids: Union[UUID, List[UUID]]) -> "VariantIntervalCollection":
        if isinstance(id_or_ids, UUID):
            ids = [id_or_ids]
        else:
            ids = id_or_ids
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
        """Build a :class:`VariantIntervalCollection` from a dictionary representation"""
        return VariantIntervalCollection(
            variant_intervals=[
                VariantInterval.from_dict(x, parent_or_seq_chunk_parent) for x in vals["variant_intervals"]
            ],
            variant_collection_name=vals["variant_collection_name"],
            variant_collection_id=vals["variant_collection_id"],
            qualifiers=vals["qualifiers"],
            sequence_name=vals["sequence_name"],
            sequence_guid=vals["sequence_guid"],
            guid=vals["variant_collection_guid"],
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    def to_gff(self, chromosome_relative_coordinates: bool = True) -> Iterator[GFFRow]:
        raise NotImplementedError("Cannot export Variants to GFF")

    def to_gtf(self, chromosome_relative_coordinates: bool = True) -> Iterator[GFFRow]:
        raise NotImplementedError("Cannot export Variants to GTF")

    @property
    def id(self) -> str:
        return self.variant_collection_id

    @property
    def name(self) -> str:
        return self.variant_collection_name

    @property
    def alternative_genomic_sequence(self) -> Sequence:
        """Edited version of the original sequence"""
        if not self.has_sequence:
            raise NullSequenceException("This VariantInterval has no sequence information")
        if not self._alternative_genomic_sequence:
            original_sequence = self.chunk_relative_location.parent.sequence
            original_sequence_str = str(original_sequence)

            alternative_seq_data = original_sequence_str[: self.variant_intervals[0].chunk_relative_location.start]
            for i in range(len(self.variant_intervals) - 1):
                curr_edit = self.variant_intervals[i]
                next_edit = self.variant_intervals[i + 1]
                alternative_seq_data += str(curr_edit.sequence)
                alternative_seq_data += original_sequence_str[
                    curr_edit.chunk_relative_location.end : next_edit.chunk_relative_location.start
                ]
            alternative_seq_data += str(self.variant_intervals[-1].sequence)
            alternative_seq_data += original_sequence_str[self.variant_intervals[-1].chunk_relative_location.end :]

            self._alternative_genomic_sequence = Sequence(
                data=alternative_seq_data,
                alphabet=original_sequence.alphabet,
                id=None,
                type=",".join(sorted(self.variant_types)),
                parent=None,
                validate_alphabet=False,
            )
        return self._alternative_genomic_sequence

    @property
    def parent_with_alternative_sequence(self) -> Parent:
        if not self.has_sequence:
            raise NullSequenceException("This VariantInterval has no sequence information")
        if self._parent_with_alternative_sequence is None:
            # have to import here to avoid circular imports
            from biocantor.io.parser import seq_chunk_to_parent, seq_to_parent

            if self.chunk_relative_location.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
                self._parent_with_alternative_sequence = seq_chunk_to_parent(
                    str(self.alternative_genomic_sequence),
                    self.chunk_relative_location.parent_id,
                    self.chunk_relative_location.parent.parent.location.start,
                    self.chunk_relative_location.parent.parent.location.start + len(self.alternative_genomic_sequence),
                )
            else:
                self._parent_with_alternative_sequence = seq_to_parent(str(self.alternative_genomic_sequence))
        return self._parent_with_alternative_sequence

    def lift_over_location(self, location: Location) -> Location:
        """
        Construct a new Location that takes the alternative sequence defined by this VariantIntervalCollection into
        account.

        The Location can be chunk-relative or chromosome-relative. It will be returned relative to the coordinate
        system of this VariantInterval.
        """
        if location is EmptyLocation():
            return EmptyLocation()

        # all liftovers happen in chromosome coordinates
        if location.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            location = location.lift_over_to_first_ancestor_of_type(SequenceType.CHROMOSOME)

        # lift over sequenceless to avoid overhead
        if isinstance(location, SingleInterval):
            for variant in self.variant_intervals:
                location = variant._lift_over_chromosome_location_single_interval(location)
        elif isinstance(location, CompoundInterval):
            for variant in self.variant_intervals:
                location = variant._lift_over_chromosome_location_compound_interval(location)
        else:
            raise ValueError("Invalid Location type passed")
        if self.has_sequence:
            return self.liftover_location_to_seq_chunk_parent(location, self.parent_with_alternative_sequence)
        else:
            return location
