import pytest

from biocantor.exc import EmptyLocationException, MismatchedParentException
from biocantor.location.location_impl import (
    EmptyLocation,
    SingleInterval,
    CompoundInterval,
)
from biocantor.location.strand import Strand


class TestEmptyLocation:
    def test_is_singleton(self):
        assert EmptyLocation() is EmptyLocation()

    def test_len(self):
        assert len(EmptyLocation()) == 0

    def test_length(self):
        assert EmptyLocation().length == 0

    def test_parent(self):
        assert EmptyLocation().parent is None

    def test_strand(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().strand

    def test_start(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().start

    def test_end(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().end

    def test_is_contiguous(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().is_contiguous

    def test_blocks(self):
        assert EmptyLocation().blocks == []

    def test_scan_blocks(self):
        assert not EmptyLocation().scan_blocks()

    def test_num_blocks(self):
        assert EmptyLocation().num_blocks == 0

    def test_optimize_blocks(self):
        assert EmptyLocation().optimize_blocks() == EmptyLocation()

    def test_gap_list(self):
        assert EmptyLocation().gap_list() == []

    def test_gaps_location(self):
        assert EmptyLocation().gaps_location() == EmptyLocation()

    def test_extract_sequence(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().extract_sequence()

    def test_parent_to_relative_pos(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().parent_to_relative_pos(0)

    def test_relative_to_parent_pos(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().relative_to_parent_pos(0)

    def test_parent_to_relative_location(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().parent_to_relative_location(SingleInterval(0, 0, Strand.PLUS))

    @pytest.mark.parametrize(
        "other",
        [
            EmptyLocation(),
            SingleInterval(5, 10, Strand.PLUS),
            CompoundInterval([5], [10], Strand.PLUS),
        ],
    )
    def test_intersection(self, other):
        assert EmptyLocation().intersection(other) == EmptyLocation()

    def test_intersection_error(self):
        with pytest.raises(MismatchedParentException):
            EmptyLocation().intersection(SingleInterval(0, 1, Strand.PLUS, parent="seq"), strict_parent_compare=True)

    @pytest.mark.parametrize(
        "other",
        [
            EmptyLocation(),
            SingleInterval(5, 10, Strand.PLUS),
            CompoundInterval([5], [10], Strand.PLUS),
        ],
    )
    def test_minus(self, other):
        assert EmptyLocation().minus(other) == EmptyLocation()

    def test_minus_error(self):
        with pytest.raises(MismatchedParentException):
            EmptyLocation().minus(SingleInterval(0, 1, Strand.PLUS, parent="seq"), strict_parent_compare=True)

    @pytest.mark.parametrize(
        "other",
        [
            EmptyLocation(),
            SingleInterval(5, 10, Strand.PLUS),
            CompoundInterval([5], [10], Strand.PLUS),
        ],
    )
    def test_location_relative_to(self, other):
        assert EmptyLocation().location_relative_to(other) == EmptyLocation()

    def test_relative_interval_to_parent_location(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().relative_interval_to_parent_location(0, 1, Strand.PLUS)

    def test_has_overlap(self):
        assert EmptyLocation().has_overlap(EmptyLocation()) is False

    def test_has_overlap_error(self):
        with pytest.raises(MismatchedParentException):
            EmptyLocation().has_overlap(SingleInterval(0, 1, Strand.PLUS, parent="seq"), strict_parent_compare=True)

    def test_contains(self):
        assert not EmptyLocation().contains(EmptyLocation())
        assert not EmptyLocation().contains(SingleInterval(0, 1, Strand.PLUS))

    def test_contains_error(self):
        with pytest.raises(MismatchedParentException):
            EmptyLocation().contains(SingleInterval(0, 1, Strand.PLUS, parent="seq"), strict_parent_compare=True)

    def test_reverse(self):
        assert EmptyLocation().reverse() == EmptyLocation()

    def test_reverse_strand(self):
        assert EmptyLocation().reverse_strand() == EmptyLocation()

    def test_reset_strand(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().reset_strand(Strand.PLUS)

    def test_reset_parent(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().reset_parent(None)

    def test_extend_absolute(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().extend_absolute(0, 0)

    def test_extend_relative(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().extend_relative(0, 0)

    def test_shift_position(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().shift_position(0)

    def test_distance_to(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().distance_to(EmptyLocation())

    def test_merge_overlapping(self):
        assert EmptyLocation().merge_overlapping() == EmptyLocation()

    def test_first_ancestor_of_type(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().first_ancestor_of_type("seqtype")

    def test_has_ancestor_of_type(self):
        assert EmptyLocation().has_ancestor_of_type("seqtype") is False

    def test_union(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().union(SingleInterval(0, 1, Strand.PLUS))

    def test_union_preserve_overlaps(self):
        with pytest.raises(EmptyLocationException):
            EmptyLocation().union_preserve_overlaps(SingleInterval(0, 1, Strand.PLUS))
