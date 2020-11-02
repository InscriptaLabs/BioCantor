from functools import total_ordering, reduce
from typing import Optional, List, Iterator

from Bio.SeqFeature import FeatureLocation, CompoundLocation
from inscripta.biocantor.util.types import ParentInputType
from inscripta.biocantor.exc import (
    InvalidStrandException,
    InvalidPositionException,
    UnsupportedOperationException,
    EmptyLocationException,
)
from inscripta.biocantor.location.distance import DistanceType
from inscripta.biocantor.location.location import Location
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import Parent, make_parent
from inscripta.biocantor.sequence import Sequence
from inscripta.biocantor.util.object_validation import ObjectValidation


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
        parent_obj = make_parent(parent) if parent else None
        if parent_obj and parent_obj.sequence:
            if not start < len(parent_obj.sequence):
                raise InvalidPositionException(
                    f"Start position ({start}) must be < parent length ({len(parent_obj.sequence)})"
                )
            if not end <= len(parent_obj.sequence):
                raise InvalidPositionException(
                    f"End position ({end}) must be <= parent length ({len(parent_obj.sequence)})"
                )

        self._start = start
        self._end = end
        self._strand = strand
        self._parent = parent_obj.reset_location(SingleInterval(start, end, strand)) if parent_obj else None
        self._sequence = None

    def __str__(self):
        return f"{self.start}-{self.end}:{self.strand}"

    def __repr__(self):
        data = f"{repr(self.parent)}:{str(self)}" if self.parent else str(self)
        return f"<SingleInterval {data}>"

    def __len__(self):
        return self._end - self._start

    @property
    def parent(self) -> Parent:
        return self._parent

    @property
    def strand(self) -> Strand:
        return self._strand

    @property
    def start(self) -> int:
        return self._start

    @property
    def end(self) -> int:
        return self._end

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

    def __eq__(self, other):
        if type(other) is not SingleInterval:
            return False
        return (
            self.start == other.start
            and self.end == other.end
            and self.strand == other.strand
            and self.parent == other.parent
        )

    def optimize_blocks(self) -> Location:
        if len(self) == 0:
            return EmptyLocation()
        return self

    def gap_list(self) -> List["Location"]:
        return []

    def gaps_location(self) -> "Location":
        return EmptyLocation()

    def __hash__(self):
        return hash((self.start, self.end, self.strand, self.parent))

    def __lt__(self, other: Location):
        self_parent_id = self.parent.id if self.parent is not None else ""
        other_parent_id = other.parent.id if other.parent is not None else ""
        if self_parent_id != other_parent_id:
            return self_parent_id < other_parent_id
        if self.start != other.start:
            return self.start < other.start
        if self.end != other.end:
            return self.end < other.end
        return self.strand < other.strand

    def extract_sequence(self) -> Sequence:
        if self._sequence is None:
            ObjectValidation.require_location_has_parent_with_sequence(self)
            seq_plus_strand = str(self.parent.sequence)[self.start : self.end]
            if self.strand == Strand.PLUS:
                self._sequence = Sequence(
                    seq_plus_strand,
                    self.parent.sequence.alphabet,
                    validate_alphabet=False,
                )
            elif self.strand == Strand.MINUS:
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

    def has_overlap(self, other: Location, match_strand: bool = False) -> bool:
        # Todo: Should this always be 'False'?  If our intervals are on different chromosomes, they can't overlap
        #   (chr1:100-200 doesn't overlap chr2:150-250) but just because they have different parents doesn't mean
        #   they can't overlap...
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
            return other.has_overlap(self, match_strand)

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

    def intersection(self, other: Location, match_strand: bool = True) -> Location:
        if not self.has_overlap(other, match_strand=match_strand):
            return EmptyLocation()
        if type(other) is SingleInterval:
            return self._intersection_single_interval(other)
        intersect_other_strand = other.intersection(self, match_strand=match_strand)
        return intersect_other_strand.reset_strand(self.strand)

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
        else:
            return CompoundInterval.from_single_intervals([self, other])

    def minus(self, other: Location, match_strand: bool = True) -> Location:
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

    def _location_relative_to(self, other: Location) -> Location:
        intersection = other.intersection(self, match_strand=False)
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

    def __init__(
        self,
        starts: List[int],
        ends: List[int],
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
            raise ValueError("Lists of start end end positions must be nonempty and have same length")
        starts_sorted = sorted(starts)
        ends_sorted = sorted(ends)
        parent_obj = make_parent(parent) if parent else None
        self._parent = parent_obj.reset_location(CompoundInterval(starts, ends, strand)) if parent_obj else None
        self._strand = strand
        single_interval_parent = self.parent.strip_location_info() if self.parent else None
        single_intervals_no_parent = [
            SingleInterval(starts_sorted[i], ends_sorted[i], self.strand) for i in range(len(starts_sorted))
        ]
        self._single_intervals = [
            interval.reset_parent(single_interval_parent) for interval in single_intervals_no_parent
        ]
        self._start = self._single_intervals[0].start
        self._end = self._single_intervals[-1].end

    @classmethod
    def from_single_intervals(cls, intervals: List[SingleInterval]) -> "CompoundInterval":
        errors = []
        if not intervals:
            errors.append("List of intervals must be nonempty")
        if len(set([interval.strand for interval in intervals])) > 1:
            errors.append(f"Intervals must all have same strand: {set([interval.strand for interval in intervals])}")
        interval_parents = set(
            [interval.parent.strip_location_info() if interval.parent else None for interval in intervals]
        )
        if len(interval_parents) > 1:
            errors.append(f"Intervals must all have same parent: {set([interval.parent.id for interval in intervals])}")
        if errors:
            raise ValueError("\n".join(errors))
        interval_starts = [interval.start for interval in intervals]
        interval_ends = [interval.end for interval in intervals]
        strand = intervals[0].strand
        new_parent = interval_parents.pop() if interval_parents else None
        return cls(interval_starts, interval_ends, strand, new_parent)

    def __str__(self):
        return ", ".join([str(interval) for interval in self._single_intervals])

    def __eq__(self, other):
        if type(other) is not CompoundInterval:
            return False
        return self.blocks == other.blocks

    def __repr__(self):
        return f"CompoundInterval <{str(self)}>"

    def __len__(self):
        return sum([len(interval) for interval in self._single_intervals])

    @property
    def parent(self) -> Parent:
        return self._parent

    @property
    def strand(self) -> Strand:
        return self._strand

    @property
    def start(self) -> int:
        return self._start

    @property
    def end(self) -> int:
        return self._end

    @property
    def num_blocks(self):
        return len(self._single_intervals)

    @property
    def is_contiguous(self) -> bool:
        return all(
            [self._single_intervals[i + 1].start == self._single_intervals[i].end for i in range(self.num_blocks - 1)]
        )

    @property
    def is_overlapping(self) -> bool:
        """Does this interval have overlaps?"""
        return any(
            [self._single_intervals[i].end > self._single_intervals[i + 1].start for i in range(self.num_blocks - 1)]
        )

    @property
    def is_empty(self) -> bool:
        return self == EmptyLocation()

    @property
    def blocks(self) -> List[Location]:
        return self._single_intervals

    def scan_blocks(self) -> Iterator[Location]:
        self.strand.assert_directional()
        if self.strand == Strand.PLUS:
            for block in self.blocks:
                yield block
        if self.strand == Strand.MINUS:
            for i in range(self.num_blocks):
                yield self.blocks[self.num_blocks - i - 1]

    def extract_sequence(self) -> Sequence:
        self.strand.assert_directional()
        block_seqs = [interval.extract_sequence() for interval in self._single_intervals]
        if self.strand == Strand.PLUS:
            return reduce(lambda seq1, seq2: seq1.append(seq2), block_seqs)
        else:
            return reduce(lambda seq1, seq2: seq1.append(seq2), reversed(block_seqs))

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

        def do_work(rel_pos: int, blocks: List[SingleInterval]):
            plus_strand = self.strand == Strand.PLUS
            block0 = blocks[0] if plus_strand else blocks[-1]
            block0_len = len(block0)
            if rel_pos < block0_len:
                return block0.start + rel_pos if plus_strand else block0.end - 1 - rel_pos
            blocks_tail = blocks[1:] if plus_strand else blocks[:-1]
            return do_work(rel_pos - block0_len, blocks_tail)

        return do_work(relative_pos, self._single_intervals)

    def relative_interval_to_parent_location(
        self, relative_start: int, relative_end: int, relative_strand: Strand
    ) -> Location:
        start_on_parent = self.relative_to_parent_pos(relative_start)
        end_on_parent_inclusive = self.relative_to_parent_pos(relative_end - 1)
        parent_start = min(start_on_parent, end_on_parent_inclusive)
        parent_end = max(start_on_parent, end_on_parent_inclusive) + 1
        intersect_same_strand = self.intersection(
            SingleInterval(
                parent_start,
                parent_end,
                self.strand,
                parent=self.parent.strip_location_info() if self.parent else None,
            )
        )
        return intersect_same_strand.reset_strand(relative_strand.relative_to(self.strand))

    def has_overlap(self, other: Location, match_strand: bool = False) -> bool:
        return any([interval.has_overlap(other, match_strand) for interval in self._single_intervals])

    def optimize_blocks(self) -> Location:
        return self._combine_blocks()._remove_empty_blocks()._to_single_interval_if_one_block()

    def gap_list(self) -> List["Location"]:
        optimized = self.optimize_blocks()
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

    def _remove_empty_blocks(self) -> "CompoundInterval":
        return CompoundInterval.from_single_intervals([block for block in self._single_intervals if len(block) > 0])

    def _to_single_interval_if_one_block(self) -> Location:
        return self if self.num_blocks > 1 else self._single_intervals[0]

    def _combine_blocks(self) -> "CompoundInterval":
        new_blocks = []
        curr_block = self._single_intervals[0]
        i = 1
        while i < self.num_blocks:
            next_block = self._single_intervals[i]
            if curr_block.end >= next_block.start:
                new_parent = curr_block.parent.strip_location_info() if curr_block.parent else None
                curr_block = SingleInterval(
                    curr_block.start,
                    next_block.end,
                    curr_block.strand,
                    parent=new_parent,
                )
            else:
                new_blocks.append(curr_block)
                curr_block = next_block
            i += 1
        new_blocks.append(curr_block)
        return CompoundInterval.from_single_intervals(new_blocks)

    def reverse(self) -> "CompoundInterval":
        def reflect_position(relative_pos):
            return self.start + self.end - relative_pos

        def reflect_interval(single_interval: SingleInterval) -> SingleInterval:
            start = reflect_position(single_interval.end)
            end = reflect_position(single_interval.start)
            strand = single_interval.strand.reverse()
            return SingleInterval(start, end, strand)

        new_parent = self.parent.strip_location_info() if self.parent else None
        return CompoundInterval.from_single_intervals(
            [reflect_interval(single_interval) for single_interval in self._single_intervals]
        ).reset_parent(new_parent)

    def reverse_strand(self) -> "CompoundInterval":
        return CompoundInterval.from_single_intervals(
            [interval.reverse_strand() for interval in self._single_intervals]
        )

    def reset_strand(self, new_strand: Strand) -> Location:
        return CompoundInterval.from_single_intervals(
            [interval.reset_strand(new_strand) for interval in self._single_intervals]
        )

    def reset_parent(self, new_parent: Parent) -> "CompoundInterval":
        return CompoundInterval.from_single_intervals(
            [interval.reset_parent(new_parent) for interval in self._single_intervals]
        )

    def shift_position(self, shift: int) -> Location:
        return CompoundInterval.from_single_intervals(
            [interval.shift_position(shift) for interval in self._single_intervals]
        )

    def distance_to(self, other: Location, distance_type: DistanceType = DistanceType.INNER) -> int:
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

    def intersection(self, other: Location, match_strand: bool = True) -> Location:
        if not self.has_overlap(other, match_strand=match_strand):
            return EmptyLocation()
        if type(other) is SingleInterval:
            return self._intersection_single_interval(other, match_strand=match_strand)
        if type(other) is CompoundInterval:
            return self._intersection_compound_interval(other, match_strand=match_strand)
        raise UnsupportedOperationException(f"Not implemented for type {type(other)}")

    def _intersection_single_interval(self, other: Location, match_strand: bool) -> Location:
        ObjectValidation.require_object_has_type(other, SingleInterval)
        interval_intersections = []
        for single_interval in self._single_intervals:
            if single_interval.has_overlap(other, match_strand=match_strand):
                interval_intersections.append(single_interval.intersection(other, match_strand=match_strand))
        return CompoundInterval.from_single_intervals(interval_intersections).optimize_blocks()

    def _intersection_compound_interval(self, other: Location, match_strand: bool) -> Location:
        ObjectValidation.require_object_has_type(other, CompoundInterval)
        interval_intersections = []
        for self_single_interval in self._single_intervals:
            if self_single_interval.has_overlap(other, match_strand=match_strand):
                for other_single_interval in other.blocks:
                    if self_single_interval.has_overlap(other_single_interval, match_strand=match_strand):
                        interval_intersections.append(
                            self_single_interval.intersection(other_single_interval, match_strand=match_strand)
                        )
        return CompoundInterval.from_single_intervals(interval_intersections).optimize_blocks()

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
            return CompoundInterval.from_single_intervals(
                non_overlapping_blocks + [interval_all_overlappers]
            ).optimize_blocks()
        return CompoundInterval.from_single_intervals(non_overlapping_blocks + [other]).optimize_blocks()

    @staticmethod
    def _merge_compound_blocks(blocks: List[SingleInterval]) -> Location:
        return reduce(lambda left, right: left.union(right), blocks)

    def _union_compound_interval(self, other: Location):
        ObjectValidation.require_object_has_type(other, CompoundInterval)
        blocks = sorted(self.blocks + other.blocks)
        return self._merge_compound_blocks(blocks)

    def minus(self, other: Location, match_strand: bool = True) -> Location:
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

    def _location_relative_to(self, other: Location) -> Location:
        rel_loc_each_single_interval = [
            other.parent_to_relative_location(block) for block in self.blocks if other.has_overlap(block)
        ]
        return reduce(
            lambda location1, location2: location1.union(location2),
            rel_loc_each_single_interval,
        ).optimize_blocks()

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

    def __len__(self):
        return 0

    def __str__(self):
        return "EmptyLocation"

    def __eq__(self, other):
        return type(other) is _EmptyLocation

    def __repr__(self):
        return "EmptyLocation"

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

    def has_overlap(self, other: Location, match_strand: bool = False) -> bool:
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

    def location_relative_to(self, other: Location) -> Location:
        return self

    def _location_relative_to(self, other: Location) -> Location:
        return self

    def distance_to(self, other: Location, distance_type: DistanceType = DistanceType.INNER) -> int:
        raise EmptyLocationException

    def intersection(self, other: Location, match_strand: bool = True) -> Location:
        return self

    def union(self, other: Location) -> Location:
        raise EmptyLocationException

    def minus(self, other: Location, match_strand: bool = True) -> Location:
        return self

    def extend_absolute(self, extend_start: int, extend_end: int) -> Location:
        raise EmptyLocationException

    def extend_relative(self, extend_upstream: int, extend_downstream: int) -> Location:
        raise EmptyLocationException

    def merge_overlapping(self) -> Location:
        return self

    def to_biopython(self):
        raise EmptyLocationException

    def first_ancestor_of_type(self, sequence_type: str) -> Parent:
        raise EmptyLocationException

    _instance = None


def EmptyLocation():
    """Returns the single EmptyLocation instance"""
    if _EmptyLocation._instance is None:
        _EmptyLocation._instance = _EmptyLocation()
    return _EmptyLocation._instance
