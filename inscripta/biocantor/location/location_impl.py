from functools import total_ordering, reduce
from itertools import chain, islice
from typing import Optional, List, Iterator, Union, Tuple

from Bio.SeqFeature import FeatureLocation, CompoundLocation

from inscripta.biocantor import DistanceType
from inscripta.biocantor.exc import (
    InvalidStrandException,
    InvalidPositionException,
    UnsupportedOperationException,
    EmptyLocationException,
    LocationException,
)
from inscripta.biocantor.location.location import Location
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import Parent, make_parent, SequenceType
from inscripta.biocantor.sequence import Sequence
from inscripta.biocantor.util.object_validation import ObjectValidation
from inscripta.biocantor.util.ordering import RelativeOrder
from inscripta.biocantor.util.types import ParentInputType

try:
    import cgranges
except ImportError:
    HAS_CGRANGES = False
else:
    HAS_CGRANGES = True


@total_ordering
class SingleInterval(Location):
    """A single contiguous interval within a sequence"""

    def __init__(
        self,
        start: int,
        end: int,
        strand: Strand,
        parent: Optional[ParentInputType] = None,
    ):
        """
        Parameters
        ----------
        start
            0-based start position of this Location on parent
        end
            0-based exclusive end position of this Location on parent
        strand
            Strand of this Location on parent
        parent
            Parent containing this Location. If parent is not provided including a Sequence,
            some methods requiring the sequence cannot be called. Additionally, some position validation
            cannot be performed. If parent has location attribute, location is ignored and set to this location.
        """
        if not 0 <= start <= end:
            raise InvalidPositionException(f"Positions must satisfy 0 <= start <= end. Start: {start}, end: {end}")

        self.start = start
        self.end = end
        self.strand = strand
        self.length = end - start
        self._sequence = None
        self.parent = None

        if parent:
            parent_obj = make_parent(parent)
            if parent_obj.sequence is not None and end > len(parent_obj.sequence):
                raise InvalidPositionException(
                    f"End position ({end}) must be <= parent length ({len(parent_obj.sequence)})"
                )
            self.parent = parent_obj.reset_location(SingleInterval(start, end, strand))

    def __str__(self):
        return f"{self.start}-{self.end}:{self.strand}"

    def __repr__(self):
        data = f"{repr(self.parent)}:{str(self)}" if self.parent else str(self)
        return f"<SingleInterval {data}>"

    @property
    def is_contiguous(self) -> bool:
        return True

    @property
    def is_empty(self) -> bool:
        return self == EmptyLocation()

    @property
    def blocks(self) -> List[Location]:
        return [self]

    def scan_blocks(self) -> Iterator[Location]:
        yield self

    @property
    def num_blocks(self) -> int:
        return 1

    @property
    def is_overlapping(self) -> bool:
        """SingleInterval is by definition always non-overlapping"""
        return False

    @property
    def _full_span_interval(self) -> Location:
        return self

    def __eq__(self, other):
        if type(other) is not SingleInterval:
            return False
        if self.start != other.start:
            return False
        if self.end != other.end:
            return False
        if self.strand is not other.strand:
            return False
        if self.parent != other.parent:
            return False
        return True

    def optimize_blocks(self) -> Location:
        if len(self) == 0:
            return EmptyLocation()
        return self

    def gap_list(self) -> List["Location"]:
        return []

    def gaps_location(self) -> "Location":
        return EmptyLocation()

    def __hash__(self):
        return hash((self.start, self.end, self.strand, self.parent.id if self.parent else 0))

    def __lt__(self, other: Location):
        return self.compare(other) < 0

    def compare(self, other: Location) -> RelativeOrder:
        """Returns a negative integer if this Location is less than the other Location, positive integer if it is
        greater, and zero otherwise."""
        self_parent_id = self.parent.id if self.parent is not None else ""
        other_parent_id = other.parent.id if other.parent is not None else ""
        if self_parent_id != other_parent_id:
            return RelativeOrder.LESS if self_parent_id < other_parent_id else RelativeOrder.GREATER
        if self.start != other.start:
            return RelativeOrder.LESS if self.start < other.start else RelativeOrder.GREATER
        if self.end != other.end:
            return RelativeOrder.LESS if self.end < other.end else RelativeOrder.GREATER
        if self.strand != other.strand:
            return RelativeOrder.LESS if self.strand < other.strand else RelativeOrder.GREATER
        return RelativeOrder.NEITHER

    def extract_sequence(self) -> Sequence:
        if self._sequence is None:
            ObjectValidation.require_location_has_parent_with_sequence(self)
            seq_plus_strand = str(self.parent.sequence)[self.start : self.end]
            if self.strand is Strand.PLUS:
                self._sequence = Sequence(
                    seq_plus_strand,
                    self.parent.sequence.alphabet,
                    validate_alphabet=False,
                )
            elif self.strand is Strand.MINUS:
                self._sequence = Sequence(
                    seq_plus_strand,
                    self.parent.sequence.alphabet,
                    validate_alphabet=False,
                ).reverse_complement()
            else:
                raise InvalidStrandException(f"Invalid strand: {self.strand}")
        return self._sequence

    def parent_to_relative_pos(self, parent_pos: int) -> int:
        if parent_pos < self.start or parent_pos >= self.end:
            raise InvalidPositionException(
                f"Requested parent position {parent_pos} lies outside this location ({self.start}-{self.end})"
            )
        if self.strand == Strand.PLUS:
            return parent_pos - self.start
        if self.strand == Strand.MINUS:
            return self.end - parent_pos - 1
        raise InvalidStrandException(f"Invalid strand: {self.strand}")

    def relative_to_parent_pos(self, relative_pos: int) -> int:
        if relative_pos < 0 or relative_pos >= len(self):
            raise ValueError(
                f"Requested position ({relative_pos}) is not between 0 and length of this location ({len(self)})"
            )
        if self.strand == Strand.PLUS:
            return self.start + relative_pos
        if self.strand == Strand.MINUS:
            return self.end - relative_pos - 1
        raise InvalidStrandException(f"Invalid strand: {self.strand}")

    def relative_interval_to_parent_location(
        self, relative_start: int, relative_end: int, relative_strand: Strand
    ) -> Location:
        if not 0 <= relative_start <= relative_end <= len(self):
            raise ValueError(
                f"Invalid relative interval for interval of length {len(self)}: {relative_start}-{relative_end}"
            )
        if self.strand == Strand.PLUS:
            parent_start = self.start + relative_start
            parent_end = self.start + relative_end
        elif self.strand == Strand.MINUS:
            parent_start = self.end - relative_end
            parent_end = self.end - relative_start
        else:
            raise InvalidStrandException(f"Location strand must be {Strand.PLUS} or {Strand.MINUS}")
        new_parent = self.parent.strip_location_info() if self.parent else None
        return SingleInterval(
            parent_start,
            parent_end,
            self.strand.relative_to(relative_strand),
            parent=new_parent,
        )

    def has_overlap(
        self, other: Location, match_strand: bool = False, full_span: bool = False, strict_parent_compare: bool = False
    ) -> bool:
        """Compares the overlap of this interval to another interval. If ``full_span`` is ``True``,
        then this interval is compared to the full span of the other interval, regardless of type of
        the other interval.
        """
        if strict_parent_compare:
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        if self.parent_id != other.parent_id:
            return False
        if self.parent or other.parent:
            if not (self.parent and other.parent):
                return False
            if not self.parent.equals_except_location(other.parent):
                return False
        if match_strand and self.strand != other.strand:
            return False
        if type(other) is SingleInterval:
            return self._has_overlap_single_interval(other)
        else:
            return other.has_overlap(self, match_strand, full_span)

    def _has_overlap_single_interval(self, other: Location) -> bool:
        ObjectValidation.require_object_has_type(other, SingleInterval)
        if len(self) == 0 or len(other) == 0:
            return False
        if other.start <= self.start < other.end:
            return True
        if other.start < self.end <= other.end:
            return True
        if self.start <= other.start < self.end:
            return True
        if self.start < other.end <= self.end:
            return True
        return False

    def reverse(self) -> "SingleInterval":
        return self.reverse_strand()

    def reverse_strand(self) -> "SingleInterval":
        return self.reset_strand(self.strand.reverse())

    def reset_strand(self, new_strand: Strand) -> "SingleInterval":
        new_parent = self.parent.strip_location_info() if self.parent else None
        return SingleInterval(self.start, self.end, new_strand, parent=new_parent)

    def reset_parent(self, new_parent: Parent) -> "SingleInterval":
        parent = new_parent.strip_location_info() if new_parent else None
        return SingleInterval(self.start, self.end, self.strand, parent)

    def shift_position(self, shift: int) -> "SingleInterval":
        return SingleInterval(self.start + shift, self.end + shift, self.strand, self.parent)

    def distance_to(self, other: Location, distance_type: DistanceType = DistanceType.INNER) -> int:
        ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        if distance_type == DistanceType.STARTS:
            return abs(self.start - other.start)
        if distance_type == DistanceType.ENDS:
            return abs(self.end - other.end)
        if type(other) is SingleInterval:
            return self._distance_to_single_interval(other, distance_type)
        else:
            return other.distance_to(self, distance_type)

    def _distance_to_single_interval(self, other: Location, distance_type: DistanceType) -> int:
        ObjectValidation.require_object_has_type(other, SingleInterval)
        distances = [abs(self.start - other.end), abs(self.end - other.start)]
        if distance_type == DistanceType.INNER:
            if self.has_overlap(other):
                return 0
            return min(distances)
        if distance_type == DistanceType.OUTER:
            return max(distances)
        raise NotImplementedError(f"Distance type not implemented: {distance_type.value}")

    def intersection(
        self, other: Location, match_strand: bool = True, full_span: bool = False, strict_parent_compare: bool = False
    ) -> Location:
        """Intersects this SingleInterval with another Location.

        Args:
            other: The other Location.
            match_strand: Match strand or ignore strand?
            full_span: Perform comparison on the full span of the other interval? Trivial for this SingleInterval,
                but relevant if ``other`` is a CompoundInterval.

        """
        if strict_parent_compare:
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        if not self.has_overlap(other, match_strand=match_strand, full_span=full_span):
            return EmptyLocation()
        if type(other) is SingleInterval:
            return self._intersection_single_interval(other)
        intersect_other_strand = other.intersection(self, match_strand=match_strand, full_span=full_span)
        if intersect_other_strand.strand != self.strand:
            return intersect_other_strand.reset_strand(self.strand)
        else:
            return intersect_other_strand

    def _intersection_single_interval(self, other: Location) -> Location:
        ObjectValidation.require_object_has_type(other, SingleInterval)
        new_start = max(self.start, other.start)
        new_end = min(self.end, other.end)
        new_parent = self.parent.strip_location_info() if self.parent else None
        return SingleInterval(new_start, new_end, self.strand, parent=new_parent)

    def union(self, other: Location) -> Location:
        if self.strand != other.strand:
            raise ValueError(f"Strands do not match: {self.strand} != {other.strand}")
        if self.parent:
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        if type(other) is SingleInterval:
            return self._union_single_interval(other)
        return other.union(self)

    def _union_single_interval(self, other: Location) -> Location:
        ObjectValidation.require_object_has_type(other, SingleInterval)
        new_parent = self.parent.strip_location_info() if self.parent else None
        if len(self) == 0:
            return SingleInterval(other.start, other.end, other.strand, new_parent)
        if len(other) == 0:
            return SingleInterval(self.start, self.end, self.strand, new_parent)
        if self.has_overlap(other):
            return SingleInterval(
                min(self.start, other.start),
                max(self.end, other.end),
                self.strand,
                new_parent,
            )
        if self.strand != other.strand:
            raise ValueError("Intervals must all have same strand")
        return CompoundInterval((self.start, other.start), (self.end, other.end), self.strand, new_parent)

    def union_preserve_overlaps(self, other: "Location") -> "Location":
        return _union_preserve_overlaps(self, other)

    def minus(self, other: Location, match_strand: bool = True, strict_parent_compare: bool = False) -> Location:
        if strict_parent_compare:
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        if not self.has_overlap(other, match_strand=match_strand):
            return self
        result_starts = []
        result_ends = []
        curr_result_block_start = self.start
        curr_result_block_end = self.end
        for block in other.blocks:
            if block.contains(self, match_strand=match_strand):
                return EmptyLocation()
            if block.end <= curr_result_block_start:
                continue
            if block.start >= curr_result_block_end:
                break
            if block.start >= curr_result_block_start and block.end <= curr_result_block_end:
                result_starts.append(curr_result_block_start)
                result_ends.append(block.start)
                curr_result_block_start = block.end
            elif block.start < curr_result_block_start and block.end <= curr_result_block_end:
                curr_result_block_start = block.end
            else:
                curr_result_block_end = block.start
                break
        result_starts.append(curr_result_block_start)
        result_ends.append(curr_result_block_end)
        new_parent = self.parent.strip_location_info() if self.parent else None
        return CompoundInterval(result_starts, result_ends, self.strand, parent=new_parent).optimize_blocks()

    def extend_absolute(self, extend_start: int, extend_end: int) -> Location:
        if min(extend_start, extend_end) < 0:
            raise ValueError("Extension distances must be non-negative")
        return SingleInterval(self.start - extend_start, self.end + extend_end, self.strand, self.parent)

    def extend_relative(self, extend_upstream: int, extend_downstream: int) -> Location:
        self.strand.assert_directional()
        if self.strand == Strand.PLUS:
            return self.extend_absolute(extend_upstream, extend_downstream)
        else:
            return self.extend_absolute(extend_downstream, extend_upstream)

    def _location_relative_to(
        self, other: Location, strict_parent_compare: bool = False, optimize_blocks: bool = True
    ) -> Location:
        """``optimize_blocks`` is not used here, but is still a keyword argument to ensure a unified
        API between SingleInterval and CompoundInterval."""
        intersection = other.intersection(self, match_strand=False, strict_parent_compare=strict_parent_compare)
        rel_pos_1 = other.parent_to_relative_pos(intersection.start)
        rel_pos_2 = other.parent_to_relative_pos(intersection.end - 1)
        rel_start, rel_end = min(rel_pos_1, rel_pos_2), max(rel_pos_1, rel_pos_2) + 1
        rel_strand = self.strand.relative_to(other.strand)
        return SingleInterval(rel_start, rel_end, rel_strand)

    def merge_overlapping(self) -> Location:
        return self

    def to_feature_location(self) -> FeatureLocation:
        """Convert to a BioPython FeatureLocation."""
        return FeatureLocation(self.start, self.end, self.strand.value)

    def to_biopython(self) -> FeatureLocation:
        """Provide a shared function signature with other Locations"""
        return self.to_feature_location()


class CompoundInterval(Location):
    """A location consisting of multiple intervals"""

    __slots__ = ["_single_interval_store", "_is_overlapping", "_starts", "_ends"]

    def __init__(
        self,
        starts: Union[Tuple[int, ...], List[int]],
        ends: Union[Tuple[int, ...], List[int]],
        strand: Strand,
        parent: Optional[ParentInputType] = None,
    ):
        """
        Parameters
        ----------
        starts
            0-based start positions of intervals; need not be sorted but order must match `ends`
        ends
            0-based exclusive end positions of intervals; need not be sorted but order must match `starts`
        strand
            Strand of this Location on parent
        parent
            Parent containing this Location. If parent is not provided including a Sequence,
            some methods requiring the sequence cannot be called. Additionally, some position validation
            cannot be performed. If parent has a location attribute, it is ignored and reset to this location.
        """
        if not len(starts) == len(ends) > 0:
            raise LocationException("Lists of start and end positions must be nonempty and have same length")
        parent_obj = make_parent(parent) if parent else None
        if parent_obj:
            if parent_obj.location:
                single_interval_parent = Parent(
                    id=parent_obj.id,
                    sequence_type=parent_obj.sequence_type,
                    sequence=parent_obj.sequence,
                    parent=parent_obj.parent,
                    location=CompoundInterval(starts, ends, strand),
                )
            else:
                single_interval_parent = parent_obj.reset_location(CompoundInterval(starts, ends, strand))
            self.parent = single_interval_parent
        else:
            self.parent = None
        self.strand = strand
        self._starts, self._ends = self._sort_starts_ends(starts, ends, strand)
        self._single_interval_store = None
        self._is_overlapping = None

        length = 0
        for start, end in zip(self._starts, self._ends):
            if start > end:
                raise InvalidPositionException("Block starts must be less than block ends")
            length += end - start
        self.length = length

        self.start = self._starts[0]
        self.end = self._ends[-1]

    @staticmethod
    def _sort_starts_ends(
        starts: List[int], ends: List[int], strand: Strand
    ) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
        """
        Given an array of positions and an orientation, sort the positions such that they are in incrementing order
        *relative to the orientation*.
        """
        blocks = ((start, end) for start, end in zip(starts, ends))
        if strand == Strand.PLUS:
            blocks = sorted(blocks, key=lambda x: (x[0], x[1]))
        else:
            blocks = sorted(blocks, key=lambda x: (x[0], -x[1]))
        return tuple(zip(*blocks))

    @property
    def _single_intervals(self):
        """Lazy evaluation; cached result"""
        if self._single_interval_store is None:
            self._single_interval_store = [
                SingleInterval(self._starts[i], self._ends[i], self.strand, self.parent) for i in range(self.num_blocks)
            ]
        return self._single_interval_store

    @classmethod
    def from_single_intervals(cls, intervals: List[SingleInterval]) -> "CompoundInterval":
        errors = []
        if not intervals:
            errors.append("List of intervals must be nonempty")
        if len({interval.strand for interval in intervals}) > 1:
            errors.append(f"Intervals must all have same strand: {set([interval.strand for interval in intervals])}")
        interval_parents = {
            interval.parent.strip_location_info() if interval.parent else None for interval in intervals
        }
        if len(interval_parents) > 1:
            errors.append(f"Intervals must all have same parent: {set([interval.parent.id for interval in intervals])}")
        if errors:
            raise ValueError("\n".join(errors))
        return cls._from_single_intervals_no_validation(intervals)

    @classmethod
    def _from_single_intervals_no_validation(cls, intervals: List[SingleInterval]) -> "CompoundInterval":
        starts = tuple(interval.start for interval in intervals)
        ends = tuple(interval.end for interval in intervals)
        strand = intervals[0].strand
        parent = intervals[0].parent.strip_location_info() if intervals[0].parent else None
        return cls(starts, ends, strand, parent)

    def __str__(self):
        return ", ".join([str(interval) for interval in self._single_intervals])

    def __eq__(self, other):
        if type(other) is not CompoundInterval:
            return False
        if self.num_blocks != other.num_blocks:
            return False
        return all(block1 == block2 for block1, block2 in zip(self.blocks, other.blocks))

    def __hash__(self):
        return hash((self._starts, self._ends, self.strand, self.parent.id if self.parent else 0))

    def __repr__(self):
        return f"CompoundInterval <{str(self)}>"

    @property
    def num_blocks(self):
        return len(self._starts)

    @property
    def is_contiguous(self) -> bool:
        return all(next_start == end for next_start, end in zip(islice(self._starts, 1, None), self._ends))

    @property
    def is_overlapping(self) -> bool:
        """Does this interval have overlaps?"""
        if self._is_overlapping is None:
            self._is_overlapping = any(
                end > next_start for next_start, end in zip(islice(self._starts, 1, None), self._ends)
            )
        return self._is_overlapping

    @property
    def is_empty(self) -> bool:
        return self == EmptyLocation()

    @property
    def blocks(self) -> List[Location]:
        return self._single_intervals

    @property
    def _full_span_interval(self) -> SingleInterval:
        return SingleInterval(self.start, self.end, self.strand, self.parent)

    def scan_blocks(self) -> Iterator[SingleInterval]:
        self.strand.assert_directional()
        if self.strand == Strand.PLUS:
            yield from self.blocks
        if self.strand == Strand.MINUS:
            yield from reversed(self.blocks)

    def extract_sequence(self) -> Sequence:
        self.strand.assert_directional()
        if self.strand == Strand.PLUS:
            block_seqs = (interval.extract_sequence() for interval in self._single_intervals)
            return reduce(lambda seq1, seq2: seq1.append(seq2), block_seqs)
        else:
            reverse_block_seqs = (interval.extract_sequence() for interval in reversed(self._single_intervals))
            return reduce(lambda seq1, seq2: seq1.append(seq2), reverse_block_seqs)

    def parent_to_relative_pos(self, parent_pos: int) -> int:
        rel_pos = 0
        for block in self.scan_blocks():
            try:
                rel_pos += block.parent_to_relative_pos(parent_pos)
                return rel_pos
            except InvalidPositionException:
                rel_pos += len(block)
        raise InvalidPositionException(
            f"Requested parent position ({parent_pos}) lies outside this location ({str(self)})"
        )

    def relative_to_parent_pos(self, relative_pos: int) -> int:
        self.strand.assert_directional()
        if not 0 <= relative_pos < len(self):
            raise InvalidPositionException(
                f"Invalid relative position {relative_pos} for location of length {len(self)}"
            )

        if self.strand == Strand.MINUS:
            starts = reversed(self._starts)
            ends = reversed(self._ends)
            plus_strand = False
        else:
            starts = iter(self._starts)
            ends = iter(self._ends)
            plus_strand = True

        def do_work(rel_pos: int, _starts: Iterator[int], _ends: Iterator[int]):
            start, end = (next(_starts), next(_ends))
            block0_len = end - start
            if rel_pos < block0_len:
                return start + rel_pos if plus_strand else end - 1 - rel_pos
            return do_work(rel_pos - block0_len, _starts, _ends)

        return do_work(relative_pos, starts, ends)

    def relative_interval_to_parent_location(
        self, relative_start: int, relative_end: int, relative_strand: Strand
    ) -> Location:
        # validate inputs
        if relative_start > relative_end:
            raise InvalidPositionException("Relative start must be <= relative end")
        elif relative_start < 0:
            raise InvalidPositionException("Relative start must be a positive value")
        elif relative_end > len(self):
            raise InvalidPositionException("Relative end must be within the size of the interval")
        # if start == end, then just return a simple interval
        elif relative_start == relative_end:
            start_on_parent = self.relative_to_parent_pos(relative_start)
            return SingleInterval(
                start_on_parent,
                start_on_parent,
                relative_strand.relative_to(self.strand),
                parent=self.parent.strip_location_info() if self.parent else None,
            )

        # the number of relative bases we need to iterate through before starting
        remaining_len_till_start = relative_start
        # the number of relative bases we have left to consume before stopping
        remaining_len_till_end = relative_end - relative_start
        new_blocks = []
        for block in self.scan_blocks():
            # this block is to the left of the relative start, so decrement remaining_len_till_start and continue
            if len(block) <= remaining_len_till_start:
                remaining_len_till_start -= len(block)
                continue
            # this block is within the window; remaining_len_till_start might be non-zero the first time
            # the loop hits this else block
            new_sub_block_start = remaining_len_till_start
            # if this block is smaller than the remaining bases, cut the value down to size
            new_sub_block_end = min(len(block), remaining_len_till_start + remaining_len_till_end)
            # lift over this block and keep it
            sub_block = block.relative_interval_to_parent_location(new_sub_block_start, new_sub_block_end, Strand.PLUS)
            new_blocks.append(sub_block)
            # all subsequent iterations will have remaining_len_till_start == 0
            remaining_len_till_start = 0
            # decrement remaining_len_till_end until nothing is left
            remaining_len_till_end -= len(sub_block)
            # remaining_len_till_end has nothing left, so exit the algorithm
            if remaining_len_till_end < 1:
                break

        # only run reset_strand() if the strands are not the same to save object construction overhead
        new_strand = relative_strand.relative_to(self.strand)
        relative_interval = CompoundInterval._from_single_intervals_no_validation(new_blocks).optimize_blocks()
        if new_strand != self.strand:
            relative_interval = relative_interval.reset_strand(new_strand)
        return relative_interval

    def has_overlap(
        self, other: Location, match_strand: bool = False, full_span: bool = False, strict_parent_compare: bool = False
    ) -> bool:
        """If full_span is ``True``, then the full span of both this location *and* the ``other`` location
        are used for the comparison.
        """
        if strict_parent_compare:
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        if full_span:
            return self._full_span_interval.has_overlap(other, match_strand, full_span=True)
        return any((interval.has_overlap(other, match_strand, full_span=False) for interval in self._single_intervals))

    def optimize_blocks(self) -> Location:
        """
        - Removes empty blocks
        - Combines adjacent blocks, preserving strictly overlapping blocks
        - Converts to SingleInterval if has only one block
        """
        combined = self._combine_blocks(preserve_overlappers=True)
        if not combined.is_empty:
            return combined._to_single_interval_if_one_block()
        else:
            return combined

    def optimize_and_combine_blocks(self) -> Location:
        """
        - Removes empty blocks
        - Combines adjacent and overlapping blocks
        - Converts to SingleInterval if has only one block
        """
        combined = self._combine_blocks(preserve_overlappers=False)
        if not combined.is_empty:
            return combined._to_single_interval_if_one_block()
        else:
            return combined

    def gap_list(self) -> List[SingleInterval]:
        optimized = self.optimize_and_combine_blocks()
        block_iter = optimized.scan_blocks()
        gaps = []
        block1 = next(block_iter)
        for block2 in block_iter:
            gap_start = min(block1.end, block2.end)
            gap_end = max(block1.start, block2.start)
            gaps.append(SingleInterval(gap_start, gap_end, self.strand, self.parent))
            block1 = block2
        return gaps

    def gaps_location(self) -> "Location":
        gap_list = self.gap_list()
        if gap_list:
            return CompoundInterval.from_single_intervals(gap_list)
        return EmptyLocation()

    def _to_single_interval_if_one_block(self) -> Location:
        return self if self.num_blocks > 1 else self._single_intervals[0]

    def _combine_blocks(self, preserve_overlappers: bool) -> "CompoundInterval":
        """Combine adjacent and (optionally) overlapping blocks. Also strips out empty blocks.

        Parameters
        ----------
        preserve_overlappers
            Do not combine strictly overlapping blocks
        """
        new_starts = []
        new_ends = []
        curr_start = curr_end = None
        # keeps track on whether there is actually any work to be done here
        needs_combining = False
        # iterate over combined iterator of starts/ends to avoid instantiating self._blocks
        for block_start, block_end in zip(self._starts, self._ends):
            # don't include empty blocks
            if block_end - block_start == 0:
                needs_combining = True
                continue
            # base case
            elif curr_start is None:
                curr_start = block_start
                curr_end = block_end
                new_starts.append(curr_start)
                new_ends.append(curr_end)
                continue
            else:
                next_start = block_start
                next_end = block_end
                combine = (curr_end == next_start) if preserve_overlappers else (curr_end >= next_start)
                if combine:
                    new_ends[-1] = next_end
                    curr_end = next_end
                    needs_combining = True
                else:
                    new_starts.append(block_start)
                    new_ends.append(block_end)
                    curr_start = block_start
                    curr_end = block_end
        # do not build a new CompoundInterval if we don't actually need to
        if not needs_combining:
            return self
        # the new_starts/new_ends list might be empty if preserve_overlappers is True, the intervals are overlapping
        # and the block combination will result in an EmptyLocation
        if not new_starts:
            return EmptyLocation()
        new_parent = self.parent.strip_location_info() if self.parent else None
        return CompoundInterval(new_starts, new_ends, self.strand, new_parent)

    def reverse(self) -> "CompoundInterval":
        def reflect_position(relative_pos):
            return self.start + self.end - relative_pos

        new_starts = [reflect_position(interval.end) for interval in self._single_intervals]
        new_ends = [reflect_position(interval.start) for interval in self._single_intervals]
        new_strand = self.strand.reverse()
        new_parent = self.parent.strip_location_info() if self.parent else None
        return CompoundInterval(new_starts, new_ends, new_strand, new_parent)

    def reverse_strand(self) -> "CompoundInterval":
        return CompoundInterval(self._starts, self._ends, self.strand.reverse(), self.parent)

    def reset_strand(self, new_strand: Strand) -> Location:
        return CompoundInterval(self._starts, self._ends, new_strand, self.parent)

    def reset_parent(self, new_parent: Parent) -> "CompoundInterval":
        return CompoundInterval(self._starts, self._ends, self.strand, new_parent)

    def shift_position(self, shift: int) -> Location:
        starts = tuple(interval.start + shift for interval in self._single_intervals)
        ends = tuple(interval.end + shift for interval in self._single_intervals)
        r = CompoundInterval(starts, ends, self.strand, self.parent)
        # force evaluation of _single_intervals to do bounds checks
        _ = r._single_intervals
        return r

    def distance_to(self, other: Location, distance_type: DistanceType = DistanceType.INNER) -> int:
        ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        if distance_type == DistanceType.STARTS:
            return abs(self.start - other.start)
        elif distance_type == DistanceType.ENDS:
            return abs(self.end - other.end)
        elif distance_type == DistanceType.OUTER:
            return max(abs(self.start - other.end), abs(self.end - other.start))
        elif distance_type == DistanceType.INNER:
            distance = None
            for interval in self._single_intervals:
                interval_distance = other.distance_to(interval)
                if distance is None or interval_distance < distance:
                    distance = interval_distance
                    if distance == 0:
                        return 0
            return distance
        else:
            raise NotImplementedError(f"Unknown distance type {distance_type.value}")

    def intersection(
        self, other: Location, match_strand: bool = True, full_span: bool = False, strict_parent_compare: bool = False
    ) -> Location:
        """Intersects this CompoundInterval with another Location.

        Args:
            other: The other Location.
            match_strand: Match strand or ignore strand?
            full_span: Perform comparison on the full span of the other interval? In all cases, the comparison
                is performed on the full span.
            strict_parent_compare: If True, parents will be compared and an exception raised if they are not equal.
                If False, mismatched parents will result in an EmptyLocation return.

        """
        if strict_parent_compare:
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        if not self.has_overlap(other, match_strand=match_strand, full_span=full_span):
            return EmptyLocation()
        if type(other) is SingleInterval:
            return self._intersection_single_interval(other, match_strand=match_strand, full_span=full_span)
        if type(other) is CompoundInterval:
            return self._intersection_compound_interval(other, match_strand=match_strand, full_span=full_span)
        raise UnsupportedOperationException(f"Not implemented for type {type(other)}")

    def _intersection_single_interval(self, other: Location, match_strand: bool, full_span: bool = False) -> Location:
        """Intersections with full span are always symmetric full span (both are considered as full span)"""
        ObjectValidation.require_object_has_type(other, SingleInterval)
        interval_intersections = []
        if full_span is False:
            for single_interval in self._single_intervals:
                if single_interval.has_overlap(other, match_strand=match_strand, full_span=False):
                    interval_intersections.append(
                        single_interval.intersection(other, match_strand=match_strand, full_span=False)
                    )
            return CompoundInterval._from_single_intervals_no_validation(interval_intersections).optimize_blocks()
        else:
            return self._full_span_interval.intersection(other, match_strand, full_span=True)

    def _intersection_compound_interval(self, other: Location, match_strand: bool, full_span: bool = False) -> Location:
        ObjectValidation.require_object_has_type(other, CompoundInterval)
        if full_span is True:
            fs = SingleInterval(self.start, self.end, self.strand, parent=self.parent)
            return fs.intersection(other, match_strand, full_span=True)
        elif match_strand is True and self.strand != other.strand:
            return EmptyLocation()
        else:
            if HAS_CGRANGES is False:
                interval_intersections = []
                for self_single_interval in self._single_intervals:
                    if self_single_interval.has_overlap(other, match_strand=match_strand, full_span=False):
                        for other_single_interval in other.blocks:
                            if self_single_interval.has_overlap(
                                other_single_interval, match_strand=match_strand, full_span=False
                            ):
                                interval_intersections.append(
                                    self_single_interval.intersection(
                                        other_single_interval, match_strand=match_strand, full_span=False
                                    )
                                )
                return CompoundInterval._from_single_intervals_no_validation(interval_intersections).optimize_blocks()
            else:
                # construct a tree from self
                tree = cgranges.cgranges()
                for i, self_single_interval in enumerate(self._single_intervals):
                    tree.add("", self_single_interval.start, self_single_interval.end, i)
                tree.index()
                interval_intersections = []
                for other_single_interval in other.blocks:
                    for _, __, idx in tree.overlap("", other_single_interval.start, other_single_interval.end):
                        interval_intersections.append(
                            self._single_intervals[idx].intersection(
                                other_single_interval, match_strand=match_strand, full_span=False
                            )
                        )
                return CompoundInterval._from_single_intervals_no_validation(interval_intersections).optimize_blocks()

    def union(self, other: Location) -> Location:
        if self.strand != other.strand:
            raise ValueError(f"Strands do not match: {self.strand} != {other.strand}")
        if self.parent:
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        if type(other) is SingleInterval:
            return self._union_single_interval(other)
        if type(other) is CompoundInterval:
            return self._union_compound_interval(other)
        raise UnsupportedOperationException(f"Not implemented for type {type(other)}")

    def _union_single_interval(self, other: Location):
        ObjectValidation.require_object_has_type(other, SingleInterval)
        overlapping_blocks = []
        non_overlapping_blocks = []
        for block in self._single_intervals:
            if block.has_overlap(other):
                overlapping_blocks.append(block)
            else:
                non_overlapping_blocks.append(block)
        if overlapping_blocks:
            start_all_overlappers = min(other.start, min([block.start for block in overlapping_blocks]))
            end_all_overlappers = max(other.end, max([block.end for block in overlapping_blocks]))
            parent_overlappers = (
                overlapping_blocks[0].parent.strip_location_info() if overlapping_blocks[0].parent else None
            )
            interval_all_overlappers = SingleInterval(
                start_all_overlappers,
                end_all_overlappers,
                self.strand,
                parent=parent_overlappers,
            )
            return CompoundInterval._from_single_intervals_no_validation(
                non_overlapping_blocks + [interval_all_overlappers]
            ).optimize_blocks()
        return CompoundInterval._from_single_intervals_no_validation(non_overlapping_blocks + [other]).optimize_blocks()

    @staticmethod
    def _merge_compound_blocks(blocks: List[SingleInterval]) -> Location:
        return reduce(lambda left, right: left.union(right), blocks)

    def _union_compound_interval(self, other: Location):
        ObjectValidation.require_object_has_type(other, CompoundInterval)
        blocks = sorted(chain(self.blocks, other.blocks))
        return self._merge_compound_blocks(blocks)

    def union_preserve_overlaps(self, other: "Location") -> "Location":
        return _union_preserve_overlaps(self, other)

    def minus(self, other: Location, match_strand: bool = True, strict_parent_compare: bool = False) -> Location:
        if strict_parent_compare:
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        if not self.has_overlap(other, match_strand=match_strand):
            return self.optimize_blocks()
        result_blocks = []
        for block in self.blocks:
            block_minus_other = block.minus(other, match_strand=match_strand)
            if block_minus_other is not EmptyLocation():
                result_blocks.extend(block_minus_other.blocks)
        if not result_blocks:
            return EmptyLocation()
        if len(result_blocks) == 1:
            return result_blocks[0]
        return CompoundInterval.from_single_intervals(result_blocks).optimize_blocks()

    def extend_absolute(self, extend_start: int, extend_end: int) -> Location:
        if min(extend_start, extend_end) < 0:
            raise ValueError("Extension distances must be non-negative")
        extended = self
        if extend_start > 0:
            left_flank = SingleInterval(self.start - extend_start, self.start, self.strand, self.parent)
            extended = self.union(left_flank)
        if extend_end > 0:
            right_flank = SingleInterval(self.end, self.end + extend_end, self.strand, self.parent)
            extended = extended.union(right_flank)
        return extended.optimize_blocks()

    def extend_relative(self, extend_upstream: int, extend_downstream: int) -> Location:
        self.strand.assert_directional()
        if self.strand == Strand.PLUS:
            return self.extend_absolute(extend_upstream, extend_downstream)
        else:
            return self.extend_absolute(extend_downstream, extend_upstream)

    def _location_relative_to(self, other: Location, optimize_blocks: bool = True) -> Location:
        rel_loc_each_single_interval = [
            other.parent_to_relative_location(block) for block in self.blocks if other.has_overlap(block)
        ]
        if not rel_loc_each_single_interval:
            return EmptyLocation()
        rel_starts = tuple(block.start for loc in rel_loc_each_single_interval for block in loc.blocks)
        rel_ends = tuple(block.end for loc in rel_loc_each_single_interval for block in loc.blocks)
        rel_strand = rel_loc_each_single_interval[0].strand
        rel_parent = rel_loc_each_single_interval[0].parent
        parent = rel_parent.strip_location_info() if rel_parent else None
        if optimize_blocks:
            return CompoundInterval(rel_starts, rel_ends, rel_strand, parent).optimize_blocks()
        else:
            return CompoundInterval(rel_starts, rel_ends, rel_strand, parent)

    def merge_overlapping(self) -> Location:
        """If this compound interval is overlapping, merge the overlaps"""
        if not self.is_overlapping:
            return self
        return self._merge_compound_blocks(self.blocks)

    def to_compound_location(self) -> CompoundLocation:
        """Convert to a BioPython CompoundLocation, or FeatureLocation if this is a 1-block interval"""
        blocks = [b.to_biopython() for b in self.blocks]
        if len(blocks) == 1:
            return blocks[0]
        return CompoundLocation(blocks)

    def to_biopython(self) -> CompoundLocation:
        """Provide a shared function signature with other Locations"""
        return self.to_compound_location()


class _EmptyLocation(Location):
    """Singleton object representing an empty location"""

    @property
    def length(self) -> int:
        return 0

    def __str__(self):
        return "EmptyLocation"

    def __eq__(self, other):
        return type(other) is _EmptyLocation

    def __repr__(self):
        return "EmptyLocation"

    def __hash__(self):
        return hash(str(self))

    @property
    def parent(self) -> Parent:
        return None

    @property
    def strand(self) -> Strand:
        raise EmptyLocationException

    @property
    def start(self) -> int:
        raise EmptyLocationException

    @property
    def end(self) -> int:
        raise EmptyLocationException

    @property
    def is_contiguous(self) -> bool:
        raise EmptyLocationException

    @property
    def is_empty(self) -> bool:
        return True

    @property
    def blocks(self) -> List[Location]:
        return []

    def scan_blocks(self) -> None:
        return

    @property
    def num_blocks(self) -> int:
        return 0

    @property
    def is_overlapping(self) -> bool:
        """EmptyLocation is by definition always non-overlapping"""
        return False

    @property
    def _full_span_interval(self) -> Location:
        return self

    def optimize_blocks(self) -> Location:
        return self

    def gap_list(self) -> List["Location"]:
        return []

    def gaps_location(self) -> "Location":
        return self

    def extract_sequence(self) -> Sequence:
        raise EmptyLocationException

    def parent_to_relative_pos(self, parent_pos: int) -> int:
        raise EmptyLocationException

    def relative_to_parent_pos(self, relative_pos: int) -> int:
        raise EmptyLocationException

    def parent_to_relative_location(self, parent_location) -> Location:
        raise EmptyLocationException

    def relative_interval_to_parent_location(
        self, relative_start: int, relative_end: int, relative_strand: Strand
    ) -> Location:
        raise EmptyLocationException

    def has_overlap(
        self, other: Location, match_strand: bool = False, full_span: bool = False, strict_parent_compare: bool = False
    ) -> bool:
        if strict_parent_compare:
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        return False

    def reverse(self) -> Location:
        return self

    def reverse_strand(self) -> Location:
        return self

    def reset_strand(self, new_strand: Strand) -> Location:
        raise EmptyLocationException

    def reset_parent(self, new_parent: Parent) -> Location:
        raise EmptyLocationException

    def shift_position(self, shift: int) -> Location:
        raise EmptyLocationException

    def location_relative_to(self, other: Location, optimize_blocks: bool = True) -> Location:
        return self

    def _location_relative_to(self, other: Location, optimize_blocks: bool = True) -> Location:
        return self

    def distance_to(self, other: Location, distance_type: DistanceType = DistanceType.INNER) -> int:
        raise EmptyLocationException

    def intersection(
        self, other: Location, match_strand: bool = True, full_span: bool = False, strict_parent_compare: bool = False
    ) -> Location:
        if strict_parent_compare:
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        return self

    def union(self, other: Location) -> Location:
        raise EmptyLocationException

    def union_preserve_overlaps(self, other: Location) -> Location:
        raise EmptyLocationException

    def minus(self, other: Location, match_strand: bool = True, strict_parent_compare: bool = False) -> Location:
        if strict_parent_compare:
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        return self

    def extend_absolute(self, extend_start: int, extend_end: int) -> Location:
        raise EmptyLocationException

    def extend_relative(self, extend_upstream: int, extend_downstream: int) -> Location:
        raise EmptyLocationException

    def merge_overlapping(self) -> Location:
        return self

    def to_biopython(self):
        raise EmptyLocationException

    def first_ancestor_of_type(self, sequence_type: Union[str, SequenceType]) -> Parent:
        raise EmptyLocationException

    _instance = None


def EmptyLocation():
    """Returns the single EmptyLocation instance"""
    if _EmptyLocation._instance is None:
        _EmptyLocation._instance = _EmptyLocation()
    return _EmptyLocation._instance


def _union_preserve_overlaps(loc1: Location, loc2: Location) -> Location:
    if loc1.strand != loc2.strand:
        raise InvalidStrandException(f"Strands do not match: {loc1.strand} != {loc2.strand}")
    if loc1.parent:
        ObjectValidation.require_parents_equal_except_location(loc1.parent, loc2.parent)
    starts = [block.start for loc in [loc1, loc2] for block in loc.blocks]
    ends = [block.end for loc in [loc1, loc2] for block in loc.blocks]
    parent = loc1.parent.strip_location_info() if loc1.parent else None
    return CompoundInterval(starts, ends, loc1.strand, parent).optimize_blocks()
