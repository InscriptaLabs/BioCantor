import pytest
from Bio.SeqFeature import FeatureLocation, ExactPosition

from inscripta.biocantor.exc import (
    InvalidStrandException,
    NoSuchAncestorException,
    InvalidPositionException,
    MismatchedParentException,
    LocationOverlapException,
    NullParentException,
)
from inscripta.biocantor import DistanceType
from inscripta.biocantor.location.location_impl import (
    SingleInterval,
    CompoundInterval,
    EmptyLocation,
)
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import Parent
from inscripta.biocantor.sequence import Sequence
from inscripta.biocantor.sequence.alphabet import Alphabet


class TestSingleInterval:
    @pytest.mark.parametrize(
        "start,end,strand,parent,expected_exception",
        [
            # Start greater than end
            (1, 0, Strand.UNSTRANDED, None, InvalidPositionException),
            # Negative coordinate
            (-1, 1, Strand.UNSTRANDED, None, InvalidPositionException),
            # Parent too short
            (
                0,
                5,
                Strand.UNSTRANDED,
                Parent(sequence=Sequence("AAA", Alphabet.NT_STRICT)),
                InvalidPositionException,
            ),
        ],
    )
    def test_init_invalid_params(self, start, end, strand, parent, expected_exception):
        with pytest.raises(expected_exception):
            SingleInterval(start, end, strand, parent)

    @pytest.mark.parametrize(
        "interval,expected",
        [
            (SingleInterval(3, 5, Strand.MINUS, None), 2),
            (SingleInterval(1, 1, Strand.PLUS, None), 0),
        ],
    )
    def test_len(self, interval, expected):
        assert len(interval) == expected

    @pytest.mark.parametrize(
        "interval,expected",
        [
            (SingleInterval(3, 5, Strand.MINUS, None), 2),
            (SingleInterval(1, 1, Strand.PLUS, None), 0),
        ],
    )
    def test_length(self, interval, expected):
        assert interval.length == expected

    @pytest.mark.parametrize(
        "interval,expected",
        [
            (SingleInterval(3, 5, Strand.MINUS, None), None),
            # Parent gets location from interval coordinates
            (
                SingleInterval(
                    1,
                    1,
                    Strand.PLUS,
                    parent=Sequence("AAA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
                Parent(
                    location=SingleInterval(1, 1, Strand.PLUS),
                    sequence=Sequence("AAA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
            ),
        ],
    )
    def test_parent(self, interval, expected):
        assert interval.parent == expected

    @pytest.mark.parametrize(
        "interval,expected",
        [
            (
                SingleInterval(
                    0,
                    1,
                    Strand.UNSTRANDED,
                    parent=Sequence("AAA", Alphabet.NT_STRICT),
                ),
                None,
            ),
            (
                SingleInterval(
                    0,
                    1,
                    Strand.UNSTRANDED,
                    parent=Parent(id="parent", sequence=Sequence("AAA", Alphabet.NT_STRICT)),
                ),
                "parent",
            ),
        ],
    )
    def test_parent_id(self, interval, expected):
        assert interval.parent_id == expected

    @pytest.mark.parametrize(
        "interval,expected",
        [
            # Parent does not have sequence
            (
                SingleInterval(
                    0,
                    1,
                    Strand.UNSTRANDED,
                    parent=Parent(sequence_type="seq_chunk", sequence=None),
                ),
                "seq_chunk",
            ),
            # Parent has sequence with type
            (
                SingleInterval(
                    0,
                    1,
                    Strand.UNSTRANDED,
                    parent=Parent(
                        sequence_type=None,
                        sequence=Sequence("AAA", Alphabet.NT_STRICT, type="seq_chunk"),
                    ),
                ),
                "seq_chunk",
            ),
        ],
    )
    def test_parent_type(self, interval, expected):
        assert interval.parent_type == expected

    def test_strand(self):
        assert SingleInterval(3, 5, Strand.MINUS, None).strand == Strand.MINUS

    def test_start(self):
        assert SingleInterval(3, 5, Strand.MINUS, None).start == 3

    def test_end(self):
        assert SingleInterval(3, 5, Strand.MINUS, None).end == 5

    @pytest.mark.parametrize(
        "interval",
        [
            SingleInterval(3, 5, Strand.MINUS, None),
            SingleInterval(1, 1, Strand.PLUS, None),  # Empty interval is contiguous
        ],
    )
    def test_is_contiguous(self, interval):
        assert interval.is_contiguous

    def test_blocks(self):
        assert SingleInterval(3, 5, Strand.MINUS, None).blocks == [SingleInterval(3, 5, Strand.MINUS, None)]
        assert SingleInterval(3, 3, Strand.MINUS, None).blocks == [SingleInterval(3, 3, Strand.MINUS, None)]

    def test_scan_blocks(self):
        assert list(SingleInterval(3, 5, Strand.MINUS, None).scan_blocks()) == [
            SingleInterval(3, 5, Strand.MINUS, None)
        ]
        assert list(SingleInterval(3, 3, Strand.MINUS, None).scan_blocks()) == [
            SingleInterval(3, 3, Strand.MINUS, None)
        ]

    def test_num_blocks(self):
        assert SingleInterval(3, 5, Strand.MINUS, None).num_blocks == 1
        assert SingleInterval(3, 3, Strand.MINUS, None).num_blocks == 1

    @pytest.mark.parametrize(
        "interval,expected",
        [
            (SingleInterval(3, 5, Strand.MINUS), SingleInterval(3, 5, Strand.MINUS)),
            (
                SingleInterval(3, 5, Strand.MINUS, parent="parent"),
                SingleInterval(3, 5, Strand.MINUS, parent="parent"),
            ),
            (SingleInterval(3, 3, Strand.MINUS, parent="parent"), EmptyLocation()),
        ],
    )
    def test_optimize_blocks(self, interval, expected):
        assert interval.optimize_blocks() == expected

    @pytest.mark.parametrize(
        "interval",
        [
            SingleInterval(3, 5, Strand.MINUS),
            SingleInterval(3, 5, Strand.MINUS, parent="parent"),
            SingleInterval(3, 3, Strand.MINUS, parent="parent"),
        ],
    )
    def test_gap_list(self, interval):
        assert interval.gap_list() == []

    @pytest.mark.parametrize(
        "interval",
        [
            SingleInterval(3, 5, Strand.MINUS),
            SingleInterval(3, 5, Strand.MINUS, parent="parent"),
            SingleInterval(3, 3, Strand.MINUS, parent="parent"),
        ],
    )
    def test_gaps_location(self, interval):
        assert interval.gaps_location() == EmptyLocation()

    @pytest.mark.parametrize(
        "interval,other,expected",
        [
            # Different type
            (SingleInterval(1, 2, Strand.UNSTRANDED, None), "x", False),
            # No parent, other fields equal
            (
                SingleInterval(1, 2, Strand.UNSTRANDED, None),
                SingleInterval(1, 2, Strand.UNSTRANDED, None),
                True,
            ),
            # Different start position
            (
                SingleInterval(1, 2, Strand.UNSTRANDED, None),
                SingleInterval(2, 2, Strand.UNSTRANDED, None),
                False,
            ),
            # Different end position
            (
                SingleInterval(1, 3, Strand.UNSTRANDED, None),
                SingleInterval(1, 2, Strand.UNSTRANDED, None),
                False,
            ),
            # Different strand
            (
                SingleInterval(1, 2, Strand.PLUS, None),
                SingleInterval(1, 2, Strand.UNSTRANDED, None),
                False,
            ),
            # One has parent, one has no parent
            (
                SingleInterval(1, 2, Strand.UNSTRANDED, None),
                SingleInterval(
                    1,
                    2,
                    Strand.UNSTRANDED,
                    parent=Sequence("AAA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
                False,
            ),
            # Has all attributes including same parent
            (
                SingleInterval(
                    1,
                    2,
                    Strand.UNSTRANDED,
                    parent=Sequence("AAA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
                SingleInterval(
                    1,
                    2,
                    Strand.UNSTRANDED,
                    parent=Sequence("AAA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
                True,
            ),
        ],
    )
    def test_equals(self, interval, other, expected):
        assert (interval == other) is expected

    def test_hash_missing_parent(self):
        _ = SingleInterval(0, 5, Strand.PLUS).__hash__()

    @pytest.mark.parametrize(
        "interval1,interval2,expected",
        [
            # Equal intervals
            (
                SingleInterval(1, 2, Strand.PLUS, None),
                SingleInterval(1, 2, Strand.PLUS, None),
                False,
            ),
            # First has lower coordinates but parent ID is lexicographically greater
            (
                SingleInterval(
                    1,
                    2,
                    Strand.PLUS,
                    parent=Sequence(
                        "AAAAAA",
                        Alphabet.NT_STRICT,
                        id="seq2",
                        validate_alphabet=False,
                    ),
                ),
                SingleInterval(
                    5,
                    6,
                    Strand.MINUS,
                    parent=Sequence(
                        "AAAAAA",
                        Alphabet.NT_STRICT,
                        id="seq1",
                        validate_alphabet=False,
                    ),
                ),
                False,
            ),
            # Same strand, first has lower start and end
            (
                SingleInterval(1, 2, Strand.PLUS, None),
                SingleInterval(5, 6, Strand.PLUS, None),
                True,
            ),
            # Same coords, different strand
            (
                SingleInterval(1, 2, Strand.PLUS, None),
                SingleInterval(1, 2, Strand.MINUS, None),
                True,
            ),
            # Same coords, different strand
            (
                SingleInterval(1, 2, Strand.MINUS, None),
                SingleInterval(1, 2, Strand.UNSTRANDED, None),
                True,
            ),
            # Same strand, different end
            (
                SingleInterval(1, 2, Strand.PLUS, None),
                SingleInterval(1, 3, Strand.PLUS, None),
                True,
            ),
            # Same end, same strand, different start, first is less
            (
                SingleInterval(1, 2, Strand.PLUS, None),
                SingleInterval(2, 2, Strand.PLUS, None),
                True,
            ),
            # Same end, same strand, different start, second is less
            (
                SingleInterval(2, 2, Strand.PLUS, None),
                SingleInterval(1, 2, Strand.PLUS, None),
                False,
            ),
        ],
    )
    def test_lt_single_interval(self, interval1, interval2, expected):
        assert (interval1 < interval2) == expected
        if interval1 != interval2:
            assert (interval1 > interval2) == (not expected)

    @pytest.mark.parametrize(
        "interval1,interval2,expected",
        [
            # Equal intervals
            (
                SingleInterval(1, 2, Strand.PLUS, None),
                SingleInterval(1, 2, Strand.PLUS, None),
                0,
            ),
            # First has lower coordinates but parent ID is lexicographically greater
            (
                SingleInterval(
                    1,
                    2,
                    Strand.PLUS,
                    parent=Sequence(
                        "AAAAAA",
                        Alphabet.NT_STRICT,
                        id="seq2",
                        validate_alphabet=False,
                    ),
                ),
                SingleInterval(
                    5,
                    6,
                    Strand.MINUS,
                    parent=Sequence(
                        "AAAAAA",
                        Alphabet.NT_STRICT,
                        id="seq1",
                        validate_alphabet=False,
                    ),
                ),
                1,
            ),
            # Same strand, first has lower start and end
            (
                SingleInterval(1, 2, Strand.PLUS, None),
                SingleInterval(5, 6, Strand.PLUS, None),
                -1,
            ),
            # Same coords, different strand
            (
                SingleInterval(1, 2, Strand.PLUS, None),
                SingleInterval(1, 2, Strand.MINUS, None),
                -1,
            ),
            # Same coords, different strand
            (
                SingleInterval(1, 2, Strand.MINUS, None),
                SingleInterval(1, 2, Strand.UNSTRANDED, None),
                -1,
            ),
            # Same strand, different end
            (
                SingleInterval(1, 2, Strand.PLUS, None),
                SingleInterval(1, 3, Strand.PLUS, None),
                -1,
            ),
            # Same end, same strand, different start, first is less
            (
                SingleInterval(1, 2, Strand.PLUS, None),
                SingleInterval(2, 2, Strand.PLUS, None),
                -1,
            ),
            # Same end, same strand, different start, second is less
            (
                SingleInterval(2, 2, Strand.PLUS, None),
                SingleInterval(1, 2, Strand.PLUS, None),
                1,
            ),
        ],
    )
    def test_compare_single_interval(self, interval1, interval2, expected):
        assert interval1.compare(interval2) == expected

    @pytest.mark.parametrize(
        "interval,expected",
        [
            # Empty interval
            (
                SingleInterval(
                    0,
                    0,
                    Strand.MINUS,
                    parent=Sequence("AAA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
                Sequence("", Alphabet.NT_STRICT, validate_alphabet=False),
            ),
            # Entire sequence, plus strand
            (
                SingleInterval(
                    0,
                    6,
                    Strand.PLUS,
                    parent=Sequence("acgggt", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
                Sequence("acgggt", Alphabet.NT_STRICT, validate_alphabet=False),
            ),
            # Entire sequence, minus strand
            (
                SingleInterval(
                    0,
                    6,
                    Strand.MINUS,
                    parent=Sequence("acccgt", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
                Sequence("acgggt", Alphabet.NT_STRICT, validate_alphabet=False),
            ),
            # Part of sequence, plus strand
            (
                SingleInterval(
                    0,
                    3,
                    Strand.PLUS,
                    parent=Sequence("acgggt", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
                Sequence("acg", Alphabet.NT_STRICT, validate_alphabet=False),
            ),
            # Part of sequence, minus strand; parent gets stripped
            (
                SingleInterval(
                    0,
                    3,
                    Strand.MINUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="chr",
                        sequence=Sequence("acgggt", Alphabet.NT_STRICT, validate_alphabet=False),
                    ),
                ),
                Sequence("cgt", Alphabet.NT_STRICT, validate_alphabet=False),
            ),
        ],
    )
    def test_extract_sequence(self, interval, expected):
        assert interval.extract_sequence() == expected

    def test_extract_sequence_error(self):
        # No parent
        with pytest.raises(NullParentException):
            SingleInterval(5, 10, Strand.PLUS, None).extract_sequence()
        # Unstranded
        with pytest.raises(InvalidStrandException):
            SingleInterval(
                0,
                1,
                Strand.UNSTRANDED,
                parent=Sequence("AAAA", Alphabet.NT_STRICT, validate_alphabet=False),
            ).extract_sequence()

    @pytest.mark.parametrize(
        "interval,parent_pos,expected",
        [
            # First position of interval, plus strand
            (SingleInterval(5, 10, Strand.PLUS, None), 5, 0),
            # First position of interval, minus strand
            (SingleInterval(5, 10, Strand.MINUS, None), 5, 4),
            # Last position of interval, plus strand
            (SingleInterval(5, 10, Strand.PLUS, None), 9, 4),
            # Last position of interval, minus strand
            (SingleInterval(5, 10, Strand.MINUS, None), 9, 0),
            # Interior position, plus strand
            (SingleInterval(5, 10, Strand.PLUS, None), 6, 1),
            # Interior position, minus strand
            (SingleInterval(5, 10, Strand.MINUS, None), 6, 3),
            # Only position of length 1 interval, plus strand
            (SingleInterval(5, 6, Strand.PLUS, None), 5, 0),
            # Only position of length 1 interval, minus strand
            (SingleInterval(5, 6, Strand.MINUS, None), 5, 0),
            # Has parent
            (
                SingleInterval(
                    5,
                    6,
                    Strand.MINUS,
                    parent=Sequence("AAAAAA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
                5,
                0,
            ),
        ],
    )
    def test_parent_to_relative_pos(self, interval, parent_pos, expected):
        assert interval.parent_to_relative_pos(parent_pos) == expected

    @pytest.mark.parametrize(
        "interval,parent_pos",
        [
            (SingleInterval(5, 5, Strand.PLUS, None), 5),
            (SingleInterval(5, 10, Strand.PLUS, None), 4),
            (SingleInterval(5, 6, Strand.PLUS, None), 6),
        ],
    )
    def test_parent_to_relative_pos_error_invalid_pos(self, interval, parent_pos):
        with pytest.raises(InvalidPositionException):
            interval.parent_to_relative_pos(parent_pos)

    def test_parent_to_relative_pos_error_invalid_strand(self):
        with pytest.raises(InvalidStrandException):
            SingleInterval(5, 10, Strand.UNSTRANDED, None).parent_to_relative_pos(7)

    @pytest.mark.parametrize(
        "interval,relative_pos,expected",
        [
            # First position of relative interval on plus strand
            (SingleInterval(5, 10, Strand.PLUS, None), 0, 5),
            # Last position of relative interval on minus strand
            (SingleInterval(5, 10, Strand.MINUS, None), 4, 5),
            # Last position of relative interval on plus strand
            (SingleInterval(5, 10, Strand.PLUS, None), 4, 9),
            # Interior position plus strand
            (SingleInterval(5, 10, Strand.PLUS, None), 1, 6),
            # Interior position minus strand
            (SingleInterval(5, 10, Strand.MINUS, None), 3, 6),
            # Length 1 interval plus strand
            (SingleInterval(5, 6, Strand.PLUS, None), 0, 5),
            # Length 1 interval minus strand
            (SingleInterval(5, 6, Strand.MINUS, None), 0, 5),
            # Has parent
            (
                SingleInterval(
                    5,
                    6,
                    Strand.MINUS,
                    parent=Sequence("AAAAAA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
                0,
                5,
            ),
        ],
    )
    def test_relative_to_parent_pos(self, interval, relative_pos, expected):
        assert interval.relative_to_parent_pos(relative_pos) == expected

    @pytest.mark.parametrize(
        "interval,relative_pos",
        [
            # Negative position
            (SingleInterval(5, 10, Strand.PLUS, None), -1),
            # Position greater than interval length
            (SingleInterval(5, 10, Strand.MINUS, None), 8),
        ],
    )
    def test_relative_to_parent_pos_error_invalid_pos(self, interval, relative_pos):
        with pytest.raises(ValueError):
            interval.relative_to_parent_pos(relative_pos)

    def test_relative_to_parent_pos_error_invalid_strand(self):
        with pytest.raises(InvalidStrandException):
            SingleInterval(5, 10, Strand.UNSTRANDED, None).relative_to_parent_pos(2)

    def test_parent_to_relative_location_empty_location(self):
        assert SingleInterval(5, 10, Strand.PLUS).parent_to_relative_location(EmptyLocation()) is EmptyLocation()

    @pytest.mark.parametrize(
        "interval,parent_location,expected",
        [
            # Parent location is exact interval, plus/unstranded strands
            (
                SingleInterval(3, 6, Strand.PLUS),
                SingleInterval(3, 6, Strand.UNSTRANDED),
                SingleInterval(0, 3, Strand.UNSTRANDED),
            ),
            # Parent location overlaps interval to the right, plus/minus strands, both have parents with most attributes
            (
                SingleInterval(
                    0,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="chr",
                        location=SingleInterval(0, 10, Strand.PLUS),
                        parent="grandparent",
                    ),
                ),
                SingleInterval(
                    3,
                    20,
                    Strand.MINUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="chr",
                        location=SingleInterval(3, 20, Strand.MINUS),
                        parent="grandparent",
                    ),
                ),
                SingleInterval(
                    3,
                    10,
                    Strand.MINUS,
                ),
            ),
            # Parent location overlaps interval to the left, plus/minus strands, neither has parent
            (
                SingleInterval(3, 20, Strand.PLUS),
                SingleInterval(0, 10, Strand.MINUS),
                SingleInterval(0, 7, Strand.MINUS),
            ),
            # Parent location completely contains interval, plus/minus strands, have parents
            (
                SingleInterval(
                    4,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="chr",
                        location=SingleInterval(4, 10, Strand.PLUS),
                    ),
                ),
                SingleInterval(
                    3,
                    20,
                    Strand.MINUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="chr",
                        location=SingleInterval(3, 20, Strand.MINUS),
                    ),
                ),
                SingleInterval(0, 6, Strand.MINUS),
            ),
            # Interval completely contains parent interval, minus/unstranded strands, have parents
            (
                SingleInterval(0, 10, Strand.MINUS, parent=Parent(parent="grandparent")),
                SingleInterval(
                    8,
                    9,
                    Strand.UNSTRANDED,
                    parent=Parent(parent="grandparent"),
                ),
                SingleInterval(1, 2, Strand.UNSTRANDED),
            ),
        ],
    )
    def test_parent_to_relative_location_single_interval(self, interval, parent_location, expected):
        assert interval.parent_to_relative_location(parent_location) == expected

    @pytest.mark.parametrize(
        "interval,parent_location,expected_exception",
        [
            # Interval is empty
            (
                SingleInterval(5, 5, Strand.PLUS),
                SingleInterval(0, 10, Strand.PLUS),
                LocationOverlapException,
            ),
            # Parent location is empty
            (
                SingleInterval(0, 10, Strand.PLUS),
                SingleInterval(5, 5, Strand.PLUS),
                LocationOverlapException,
            ),
            # Interval has no parent, parent location has a parent
            (
                SingleInterval(1, 2, Strand.PLUS),
                SingleInterval(0, 5, Strand.PLUS, parent="parent"),
                MismatchedParentException,
            ),
            # Interval has a parent, parent location doesn't
            (
                SingleInterval(1, 2, Strand.PLUS, parent="parent"),
                SingleInterval(0, 5, Strand.PLUS),
                NullParentException,
            ),
            # Parents are different
            (
                SingleInterval(1, 2, Strand.PLUS, parent="parent1"),
                SingleInterval(0, 5, Strand.PLUS, parent="parent2"),
                MismatchedParentException,
            ),
            # No overlap
            (
                SingleInterval(0, 3, Strand.PLUS),
                SingleInterval(5, 10, Strand.PLUS),
                LocationOverlapException,
            ),
        ],
    )
    def test_parent_to_relative_location_error_single_interval(self, interval, parent_location, expected_exception):
        with pytest.raises(expected_exception):
            interval.parent_to_relative_location(parent_location)

    @pytest.mark.parametrize(
        "interval,parent_location,expected",
        [
            # Parent location is exact interval, plus/minus strands
            (
                SingleInterval(3, 5, Strand.PLUS),
                CompoundInterval([3], [5], Strand.MINUS),
                SingleInterval(0, 2, Strand.MINUS),
            ),
            # One block of parent location overlaps interval to the right, both have parents, minus/unstranded
            (
                SingleInterval(
                    3,
                    10,
                    Strand.MINUS,
                    Parent(
                        id="parent",
                        strand=Strand.MINUS,
                        location=SingleInterval(3, 10, Strand.MINUS),
                    ),
                ),
                CompoundInterval(
                    [5, 20],
                    [15, 25],
                    Strand.UNSTRANDED,
                    Parent(
                        id="parent",
                        strand=Strand.UNSTRANDED,
                        location=CompoundInterval([5, 20], [15, 25], Strand.UNSTRANDED),
                    ),
                ),
                SingleInterval(0, 5, Strand.UNSTRANDED),
            ),
            # One block of parent location overlaps interval to the left, minus/plus
            (
                SingleInterval(10, 20, Strand.MINUS),
                CompoundInterval([5, 30], [12, 40], Strand.PLUS),
                SingleInterval(8, 10, Strand.MINUS),
            ),
            # Two blocks of parent location overlap interval to the right, plus/minus
            (
                SingleInterval(3, 30, Strand.PLUS),
                CompoundInterval([10, 20, 40], [12, 35, 42], Strand.MINUS),
                CompoundInterval([7, 17], [9, 27], Strand.MINUS),
            ),
            # Two blocks of parent location overlap interval to the left, plus/unstranded
            (
                SingleInterval(30, 50, Strand.PLUS),
                CompoundInterval([10, 25, 45], [20, 32, 47], Strand.UNSTRANDED),
                CompoundInterval([0, 15], [2, 17], Strand.UNSTRANDED),
            ),
            # Interval completely contains parent location, minus/plus
            (
                SingleInterval(10, 30, Strand.MINUS),
                CompoundInterval([12, 22], [14, 24], Strand.PLUS),
                CompoundInterval([6, 16], [8, 18], Strand.MINUS),
            ),
            # Interval is completely contained in one block of parent location, plus/plus
            (
                SingleInterval(3, 5, Strand.PLUS),
                CompoundInterval([0, 10], [8, 20], Strand.PLUS),
                SingleInterval(0, 2, Strand.PLUS),
            ),
        ],
    )
    def test_parent_to_relative_location_compound_interval(self, interval, parent_location, expected):
        assert interval.parent_to_relative_location(parent_location) == expected

    @pytest.mark.parametrize(
        "interval,parent_location,expected_exception",
        [
            # Interval is empty
            (
                SingleInterval(5, 5, Strand.PLUS),
                CompoundInterval([0], [10], Strand.PLUS),
                LocationOverlapException,
            ),
            # Parent location is empty
            (
                SingleInterval(0, 10, Strand.PLUS),
                CompoundInterval([5], [5], Strand.PLUS),
                LocationOverlapException,
            ),
            # Interval has no parent, parent location has a parent
            (
                SingleInterval(1, 2, Strand.PLUS),
                CompoundInterval([0], [5], Strand.PLUS, parent="parent"),
                MismatchedParentException,
            ),
            # Interval has a parent, parent location doesn't
            (
                SingleInterval(1, 2, Strand.PLUS, parent="parent"),
                CompoundInterval([0], [5], Strand.PLUS),
                NullParentException,
            ),
            # Parents are different
            (
                SingleInterval(1, 2, Strand.PLUS, parent="parent1"),
                CompoundInterval([0], [5], Strand.PLUS, parent="parent2"),
                MismatchedParentException,
            ),
            # No overlap
            (
                SingleInterval(0, 3, Strand.PLUS),
                CompoundInterval([5], [10], Strand.PLUS),
                LocationOverlapException,
            ),
        ],
    )
    def test_parent_to_relative_location_error_compound_interval(self, interval, parent_location, expected_exception):
        with pytest.raises(expected_exception):
            interval.parent_to_relative_location(parent_location)

    @pytest.mark.parametrize(
        "interval,relative_start,relative_end,relative_strand,expected",
        [
            (
                SingleInterval(10, 20, Strand.PLUS),
                5,
                6,
                Strand.PLUS,
                SingleInterval(15, 16, Strand.PLUS),
            ),
            (
                SingleInterval(10, 20, Strand.PLUS),
                5,
                6,
                Strand.MINUS,
                SingleInterval(15, 16, Strand.MINUS),
            ),
            (
                SingleInterval(5, 6, Strand.MINUS),
                0,
                1,
                Strand.PLUS,
                SingleInterval(5, 6, Strand.MINUS, None),
            ),
            (
                SingleInterval(
                    5,
                    6,
                    Strand.MINUS,
                    parent=Sequence("ACGTACA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
                0,
                1,
                Strand.MINUS,
                SingleInterval(
                    5,
                    6,
                    Strand.PLUS,
                    parent=Sequence("ACGTACA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
            ),
            (
                SingleInterval(10, 20, Strand.MINUS, None),
                1,
                2,
                Strand.PLUS,
                SingleInterval(18, 19, Strand.MINUS, None),
            ),
            (
                SingleInterval(10, 20, Strand.MINUS, None),
                1,
                2,
                Strand.MINUS,
                SingleInterval(18, 19, Strand.PLUS, None),
            ),
            (
                SingleInterval(10, 20, Strand.MINUS, None),
                0,
                10,
                Strand.PLUS,
                SingleInterval(10, 20, Strand.MINUS, None),
            ),
            (SingleInterval(10, 20, Strand.PLUS), 5, 5, Strand.PLUS, SingleInterval(15, 15, Strand.PLUS)),
        ],
    )
    def test_relative_interval_to_parent_location(
        self, interval, relative_start, relative_end, relative_strand, expected
    ):
        assert interval.relative_interval_to_parent_location(relative_start, relative_end, relative_strand) == expected

    @pytest.mark.parametrize(
        "interval,start,end,strand,expected_error",
        [
            # End is greater than interval length
            (
                SingleInterval(
                    0,
                    4,
                    Strand.PLUS,
                    SingleInterval(0, 4, Strand.PLUS),
                ),
                3,
                5,
                Strand.PLUS,
                ValueError,
            ),
            # Start > end
            (
                SingleInterval(
                    0,
                    4,
                    Strand.PLUS,
                    SingleInterval(0, 4, Strand.PLUS),
                ),
                2,
                1,
                Strand.PLUS,
                ValueError,
            ),
            # Location on parent is unstranded
            (
                SingleInterval(0, 4, Strand.UNSTRANDED),
                1,
                2,
                Strand.PLUS,
                InvalidStrandException,
            ),
        ],
    )
    def test_relative_interval_to_parent_location_error(self, interval, start, end, strand, expected_error):
        with pytest.raises(expected_error):
            interval.relative_interval_to_parent_location(start, end, strand)

    def test_location_relative_to_empty_location(self):
        assert SingleInterval(5, 10, Strand.PLUS).location_relative_to(EmptyLocation()) is EmptyLocation()

    @pytest.mark.parametrize(
        "interval1,interval2,expected",
        [
            # Other location is exact interval, plus/unstranded strands
            (
                SingleInterval(3, 6, Strand.PLUS),
                SingleInterval(3, 6, Strand.UNSTRANDED),
                SingleInterval(0, 3, Strand.UNSTRANDED),
            ),
            # Other location overlaps interval to the right, plus/minus strands, both have parents with most attributes
            (
                SingleInterval(
                    0,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="chr",
                        location=SingleInterval(0, 10, Strand.PLUS),
                        parent="grandparent",
                    ),
                ),
                SingleInterval(
                    3,
                    20,
                    Strand.MINUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="chr",
                        location=SingleInterval(3, 20, Strand.MINUS),
                        parent="grandparent",
                    ),
                ),
                SingleInterval(
                    3,
                    10,
                    Strand.MINUS,
                ),
            ),
            # Other location overlaps interval to the left, plus/minus strands, neither has parent
            (
                SingleInterval(3, 20, Strand.PLUS),
                SingleInterval(0, 10, Strand.MINUS),
                SingleInterval(0, 7, Strand.MINUS),
            ),
            # Other location completely contains interval, plus/minus strands, have parents
            (
                SingleInterval(
                    4,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="chr",
                        location=SingleInterval(4, 10, Strand.PLUS),
                    ),
                ),
                SingleInterval(
                    3,
                    20,
                    Strand.MINUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="chr",
                        location=SingleInterval(3, 20, Strand.MINUS),
                    ),
                ),
                SingleInterval(0, 6, Strand.MINUS),
            ),
            # Interval completely contains other interval, minus/unstranded strands, have parents
            (
                SingleInterval(0, 10, Strand.MINUS, parent=Parent(parent="grandparent")),
                SingleInterval(
                    8,
                    9,
                    Strand.UNSTRANDED,
                    parent=Parent(parent="grandparent"),
                ),
                SingleInterval(1, 2, Strand.UNSTRANDED),
            ),
        ],
    )
    def test_location_relative_to_single_interval(self, interval1, interval2, expected):
        assert interval2.location_relative_to(interval1) == expected

    @pytest.mark.parametrize(
        "interval,other_location,expected_exception",
        [
            # Interval is empty
            (
                SingleInterval(5, 5, Strand.PLUS),
                SingleInterval(0, 10, Strand.PLUS),
                LocationOverlapException,
            ),
            # Other location is empty
            (
                SingleInterval(0, 10, Strand.PLUS),
                SingleInterval(5, 5, Strand.PLUS),
                LocationOverlapException,
            ),
            # Interval has no parent, other location has a parent
            (
                SingleInterval(1, 2, Strand.PLUS),
                SingleInterval(0, 5, Strand.PLUS, parent="parent"),
                NullParentException,
            ),
            # Interval has a parent, other location doesn't
            (
                SingleInterval(1, 2, Strand.PLUS, parent="parent"),
                SingleInterval(0, 5, Strand.PLUS),
                MismatchedParentException,
            ),
            # Parents are different
            (
                SingleInterval(1, 2, Strand.PLUS, parent="parent1"),
                SingleInterval(0, 5, Strand.PLUS, parent="parent2"),
                MismatchedParentException,
            ),
            # No overlap
            (
                SingleInterval(0, 3, Strand.PLUS),
                SingleInterval(5, 10, Strand.PLUS),
                LocationOverlapException,
            ),
        ],
    )
    def test_location_relative_to_error_single_interval(self, interval, other_location, expected_exception):
        with pytest.raises(expected_exception):
            interval.location_relative_to(other_location)

    @pytest.mark.parametrize(
        "compound_interval,single_interval,expected",
        [
            # Both have parent
            (
                CompoundInterval([0], [10], Strand.PLUS, parent="parent"),
                SingleInterval(5, 7, Strand.MINUS, parent="parent"),
                SingleInterval(5, 7, Strand.MINUS),
            ),
            # location is unstranded
            (
                CompoundInterval([0, 10], [5, 15], Strand.MINUS),
                SingleInterval(0, 12, Strand.UNSTRANDED),
                SingleInterval(3, 10, Strand.UNSTRANDED),
            ),
            # location contains entire location; same strand
            (
                CompoundInterval([1, 10], [2, 14], Strand.MINUS),
                SingleInterval(0, 25, Strand.MINUS),
                SingleInterval(0, 5, Strand.PLUS),
            ),
            # location overlaps two blocks; opposite strands
            (
                CompoundInterval([0, 10, 20], [3, 13, 23], Strand.PLUS),
                SingleInterval(7, 22, Strand.MINUS),
                SingleInterval(3, 8, Strand.MINUS),
            ),
            # location equals location
            (
                CompoundInterval([3], [5], Strand.PLUS),
                SingleInterval(3, 5, Strand.PLUS),
                SingleInterval(0, 2, Strand.PLUS),
            ),
        ],
    )
    def test_location_relative_to_compound_interval(self, compound_interval, single_interval, expected):
        assert single_interval.location_relative_to(compound_interval) == expected

    @pytest.mark.parametrize(
        "compound_interval,single_interval,expected_exception",
        [
            # Unstranded
            (
                CompoundInterval([3], [5], Strand.UNSTRANDED),
                SingleInterval(3, 5, Strand.PLUS),
                InvalidStrandException,
            ),
            # Location is outside location
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                SingleInterval(20, 25, Strand.PLUS),
                LocationOverlapException,
            ),
            # Location is in "intron"
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                SingleInterval(7, 8, Strand.PLUS),
                LocationOverlapException,
            ),
        ],
    )
    def test_location_relative_to_error_compound_interval(self, compound_interval, single_interval, expected_exception):
        with pytest.raises(expected_exception):
            single_interval.location_relative_to(compound_interval)

    @pytest.mark.parametrize(
        "interval,window_size,step_size,start_pos,expected",
        [
            # Window size same as location length
            (
                SingleInterval(5, 10, Strand.PLUS),
                5,
                1,
                0,
                [SingleInterval(5, 10, Strand.PLUS)],
            ),
            # Window size one less than location length
            (
                SingleInterval(5, 10, Strand.PLUS),
                4,
                1,
                0,
                [SingleInterval(5, 9, Strand.PLUS), SingleInterval(6, 10, Strand.PLUS)],
            ),
            # Minus strand
            (
                SingleInterval(5, 10, Strand.MINUS),
                4,
                1,
                0,
                [
                    SingleInterval(6, 10, Strand.MINUS),
                    SingleInterval(5, 9, Strand.MINUS),
                ],
            ),
            # Non-overlapping windows fill location
            (
                SingleInterval(0, 10, Strand.PLUS),
                2,
                3,
                2,
                [
                    SingleInterval(2, 4, Strand.PLUS),
                    SingleInterval(5, 7, Strand.PLUS),
                    SingleInterval(8, 10, Strand.PLUS),
                ],
            ),
            # Non-overlapping windows don't fill location
            (
                SingleInterval(0, 12, Strand.PLUS),
                2,
                3,
                2,
                [
                    SingleInterval(2, 4, Strand.PLUS),
                    SingleInterval(5, 7, Strand.PLUS),
                    SingleInterval(8, 10, Strand.PLUS),
                ],
            ),
        ],
    )
    def test_scan_windows(self, interval, window_size, step_size, start_pos, expected):
        assert list(interval.scan_windows(window_size, step_size, start_pos)) == expected

    @pytest.mark.parametrize(
        "interval,window_size,step_size,start_pos,expected_exception",
        [
            # Window size is greater than location length
            (SingleInterval(0, 10, Strand.PLUS), 20, 1, 0, ValueError),
            # Location is unstranded
            (SingleInterval(0, 10, Strand.UNSTRANDED), 1, 1, 0, InvalidStrandException),
            # Window size is 0
            (SingleInterval(0, 10, Strand.PLUS), 0, 1, 0, ValueError),
            # Step size is 0
            (SingleInterval(0, 10, Strand.PLUS), 1, 0, 0, ValueError),
            # Start pos is negative
            (SingleInterval(0, 10, Strand.PLUS), 1, 1, -1, ValueError),
            # Start pos is greater than location length
            (SingleInterval(0, 10, Strand.PLUS), 1, 1, 20, ValueError),
            # Start pos + window size is greater than location length
            (SingleInterval(0, 10, Strand.PLUS), 5, 1, 8, ValueError),
        ],
    )
    def test_scan_windows_error(self, interval, window_size, step_size, start_pos, expected_exception):
        with pytest.raises(expected_exception):
            list(interval.scan_windows(window_size, step_size, start_pos))

    @pytest.mark.parametrize(
        "interval,other,match_strand,expected",
        [
            (
                SingleInterval(0, 1, Strand.PLUS, parent="parent"),
                SingleInterval(0, 1, Strand.PLUS, parent="parent2"),
                True,
                False,
            ),
            (
                SingleInterval(
                    0,
                    1,
                    Strand.PLUS,
                    parent=Sequence("AA", Alphabet.NT_STRICT, id="parent"),
                ),
                SingleInterval(
                    0,
                    1,
                    Strand.PLUS,
                    parent=Sequence("AA", Alphabet.NT_STRICT, id="parent2"),
                ),
                True,
                False,
            ),
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(4, 6, Strand.PLUS),
                False,
                True,
            ),
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(5, 6, Strand.PLUS),
                True,
                False,
            ),
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(4, 6, Strand.MINUS),
                False,
                True,
            ),
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(4, 6, Strand.MINUS),
                True,
                False,
            ),
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(4, 6, Strand.UNSTRANDED),
                True,
                False,
            ),
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(4, 6, Strand.PLUS),
                False,
                True,
            ),
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(5, 6, Strand.PLUS),
                True,
                False,
            ),
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(4, 6, Strand.MINUS),
                False,
                True,
            ),
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(4, 6, Strand.MINUS),
                True,
                False,
            ),
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(4, 6, Strand.UNSTRANDED),
                True,
                False,
            ),
            (
                SingleInterval(0, 1, Strand.PLUS, parent="parent"),
                CompoundInterval([0], [1], Strand.PLUS, parent="parent2"),
                True,
                False,
            ),
            (
                SingleInterval(
                    0,
                    1,
                    Strand.PLUS,
                    parent=Sequence("AA", Alphabet.NT_STRICT, id="parent"),
                ),
                CompoundInterval(
                    [0],
                    [1],
                    Strand.PLUS,
                    parent=Sequence("AA", Alphabet.NT_STRICT, id="parent2"),
                ),
                True,
                False,
            ),
            (
                SingleInterval(10, 20, Strand.PLUS),
                CompoundInterval([5, 19], [7, 20], Strand.PLUS),
                True,
                True,
            ),
            (
                SingleInterval(6, 22, Strand.PLUS),
                CompoundInterval([5, 19], [7, 20], Strand.PLUS),
                True,
                True,
            ),
            (
                SingleInterval(7, 19, Strand.PLUS),
                CompoundInterval([5, 19], [7, 20], Strand.PLUS),
                True,
                False,
            ),
            (
                SingleInterval(6, 22, Strand.PLUS),
                CompoundInterval([5, 19], [7, 20], Strand.MINUS),
                False,
                True,
            ),
            (
                SingleInterval(6, 22, Strand.PLUS),
                CompoundInterval([5, 19], [7, 20], Strand.MINUS),
                True,
                False,
            ),
            (
                SingleInterval(3, 5, Strand.PLUS),
                SingleInterval(0, 20, Strand.MINUS),
                False,
                True,
            ),
            (
                SingleInterval(3, 5, Strand.PLUS),
                SingleInterval(0, 20, Strand.MINUS),
                True,
                False,
            ),
            (
                SingleInterval(0, 20, Strand.PLUS),
                SingleInterval(3, 5, Strand.PLUS),
                True,
                True,
            ),
            (
                SingleInterval(5, 5, Strand.PLUS),
                SingleInterval(0, 10, Strand.PLUS),
                False,
                False,
            ),
            (
                SingleInterval(0, 10, Strand.PLUS),
                SingleInterval(5, 5, Strand.PLUS),
                False,
                False,
            ),
        ],
    )
    def test_has_overlap(self, interval, other, match_strand, expected):
        assert interval.has_overlap(other, match_strand) == expected

    def test_has_overlap_error(self):
        with pytest.raises(MismatchedParentException):
            SingleInterval(0, 1, Strand.PLUS, parent="seq1").has_overlap(
                SingleInterval(0, 1, Strand.PLUS, parent="seq2"), strict_parent_compare=True
            )

    @pytest.mark.parametrize(
        "interval,expected",
        [
            (
                SingleInterval(5, 10, Strand.PLUS, None),
                SingleInterval(5, 10, Strand.MINUS, None),
            ),
            (
                SingleInterval(5, 10, Strand.UNSTRANDED, None),
                SingleInterval(5, 10, Strand.UNSTRANDED, None),
            ),
            (
                SingleInterval(
                    0,
                    2,
                    Strand.PLUS,
                    parent=Sequence("AAA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
                SingleInterval(
                    0,
                    2,
                    Strand.MINUS,
                    parent=Sequence("AAA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
            ),
            (SingleInterval(5, 10, Strand.PLUS), SingleInterval(5, 10, Strand.MINUS)),
        ],
    )
    def test_reverse_strand(self, interval, expected):
        assert interval.reverse_strand() == expected

    @pytest.mark.parametrize(
        "interval,expected",
        [
            (
                SingleInterval(5, 10, Strand.PLUS, None),
                SingleInterval(5, 10, Strand.MINUS, None),
            ),
            (
                SingleInterval(5, 10, Strand.UNSTRANDED, None),
                SingleInterval(5, 10, Strand.UNSTRANDED, None),
            ),
            (
                SingleInterval(
                    0,
                    2,
                    Strand.PLUS,
                    parent=Sequence("AAA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
                SingleInterval(
                    0,
                    2,
                    Strand.MINUS,
                    parent=Sequence("AAA", Alphabet.NT_STRICT, validate_alphabet=False),
                ),
            ),
            (SingleInterval(5, 10, Strand.PLUS), SingleInterval(5, 10, Strand.MINUS)),
        ],
    )
    def test_reverse(self, interval, expected):
        assert interval.reverse() == expected

    @pytest.mark.parametrize(
        "interval,new_strand,expected",
        [
            (
                SingleInterval(5, 10, Strand.MINUS, parent=Parent(id="parent", strand=Strand.MINUS)),
                Strand.PLUS,
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        strand=Strand.PLUS,
                        location=SingleInterval(5, 10, Strand.PLUS),
                    ),
                ),
            ),
            (
                SingleInterval(
                    5,
                    10,
                    Strand.MINUS,
                    parent=Parent(
                        id="parent",
                        strand=Strand.MINUS,
                        location=SingleInterval(5, 10, Strand.MINUS),
                    ),
                ),
                Strand.PLUS,
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        strand=Strand.PLUS,
                        location=SingleInterval(5, 10, Strand.PLUS),
                    ),
                ),
            ),
            (
                SingleInterval(5, 10, Strand.PLUS),
                Strand.MINUS,
                SingleInterval(5, 10, Strand.MINUS),
            ),
        ],
    )
    def test_reset_strand(self, interval, new_strand, expected):
        assert interval.reset_strand(new_strand) == expected

    @pytest.mark.parametrize(
        "interval,new_parent,expected",
        [
            (
                SingleInterval(5, 10, Strand.PLUS),
                Parent(
                    id="parent",
                    strand=Strand.MINUS,
                    location=SingleInterval(13, 14, Strand.MINUS),
                ),
                SingleInterval(5, 10, Strand.PLUS, parent="parent"),
            ),
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="unknown",
                        strand=Strand.PLUS,
                        location=SingleInterval(5, 10, Strand.PLUS),
                        parent="grandparent",
                    ),
                ),
                Parent(
                    id="new_parent",
                    sequence_type="chr",
                    strand=Strand.PLUS,
                    location=SingleInterval(5, 10, Strand.PLUS),
                ),
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(id="new_parent", sequence_type="chr"),
                ),
            ),
            (
                SingleInterval(5, 10, Strand.PLUS, parent="parent"),
                None,
                SingleInterval(5, 10, Strand.PLUS),
            ),
        ],
    )
    def test_reset_parent(self, interval, new_parent, expected):
        assert interval.reset_parent(new_parent) == expected

    @pytest.mark.parametrize(
        "interval,shift,expected",
        [
            # Positive shift
            (SingleInterval(0, 3, Strand.PLUS), 5, SingleInterval(5, 8, Strand.PLUS)),
            # Negative shift
            (
                SingleInterval(10, 13, Strand.MINUS),
                -5,
                SingleInterval(5, 8, Strand.MINUS),
            ),
            # Has parent
            (
                SingleInterval(
                    0,
                    3,
                    Strand.UNSTRANDED,
                    parent=Parent(
                        id="parent",
                        strand=Strand.UNSTRANDED,
                        location=SingleInterval(0, 3, Strand.UNSTRANDED),
                    ),
                ),
                5,
                SingleInterval(
                    5,
                    8,
                    Strand.UNSTRANDED,
                    parent=Parent(
                        id="parent",
                        strand=Strand.UNSTRANDED,
                        location=SingleInterval(5, 8, Strand.UNSTRANDED),
                    ),
                ),
            ),
        ],
    )
    def test_shift_position(self, interval, shift, expected):
        assert interval.shift_position(shift) == expected

    @pytest.mark.parametrize(
        "interval,shift,expected_exception",
        [
            # Off beginning of sequence
            (SingleInterval(0, 3, Strand.PLUS), -1, InvalidPositionException),
            # Off end of sequence
            (
                SingleInterval(
                    0,
                    3,
                    Strand.PLUS,
                    parent=Sequence("AAA", Alphabet.NT_STRICT),
                ),
                1,
                InvalidPositionException,
            ),
        ],
    )
    def test_shift_position_error(self, interval, shift, expected_exception):
        with pytest.raises(expected_exception):
            interval.shift_position(shift)

    @pytest.mark.parametrize(
        "interval,other,distance_type,expected",
        [
            (
                SingleInterval(5, 10, Strand.PLUS),
                SingleInterval(13, 20, Strand.MINUS),
                DistanceType.INNER,
                3,
            ),
            (
                SingleInterval(13, 20, Strand.UNSTRANDED),
                SingleInterval(0, 10, Strand.PLUS),
                DistanceType.INNER,
                3,
            ),
            (
                SingleInterval(0, 4, Strand.MINUS),
                SingleInterval(4, 7, Strand.MINUS),
                DistanceType.INNER,
                0,
            ),
            (
                SingleInterval(6, 8, Strand.PLUS),
                SingleInterval(0, 6, Strand.MINUS),
                DistanceType.INNER,
                0,
            ),
            (
                SingleInterval(10, 20, Strand.UNSTRANDED),
                SingleInterval(17, 30, Strand.UNSTRANDED),
                DistanceType.INNER,
                0,
            ),
            (
                SingleInterval(17, 30, Strand.UNSTRANDED),
                SingleInterval(10, 20, Strand.UNSTRANDED),
                DistanceType.INNER,
                0,
            ),
            (
                SingleInterval(0, 20, Strand.PLUS),
                SingleInterval(3, 5, Strand.PLUS),
                DistanceType.INNER,
                0,
            ),
            (
                SingleInterval(3, 5, Strand.PLUS),
                SingleInterval(0, 20, Strand.MINUS),
                DistanceType.INNER,
                0,
            ),
            (
                SingleInterval(10, 20, Strand.MINUS),
                SingleInterval(10, 15, Strand.UNSTRANDED),
                DistanceType.INNER,
                0,
            ),
            (
                SingleInterval(10, 20, Strand.MINUS),
                SingleInterval(13, 20, Strand.UNSTRANDED),
                DistanceType.INNER,
                0,
            ),
            (
                SingleInterval(5, 10, Strand.PLUS),
                SingleInterval(13, 20, Strand.MINUS),
                DistanceType.OUTER,
                15,
            ),
            (
                SingleInterval(13, 20, Strand.UNSTRANDED),
                SingleInterval(0, 10, Strand.PLUS),
                DistanceType.OUTER,
                20,
            ),
            (
                SingleInterval(0, 4, Strand.MINUS),
                SingleInterval(4, 7, Strand.MINUS),
                DistanceType.OUTER,
                7,
            ),
            (
                SingleInterval(6, 8, Strand.PLUS),
                SingleInterval(0, 6, Strand.MINUS),
                DistanceType.OUTER,
                8,
            ),
            (
                SingleInterval(10, 20, Strand.UNSTRANDED),
                SingleInterval(17, 30, Strand.UNSTRANDED),
                DistanceType.OUTER,
                20,
            ),
            (
                SingleInterval(17, 30, Strand.UNSTRANDED),
                SingleInterval(10, 20, Strand.UNSTRANDED),
                DistanceType.OUTER,
                20,
            ),
            (
                SingleInterval(0, 20, Strand.PLUS),
                SingleInterval(3, 5, Strand.PLUS),
                DistanceType.OUTER,
                17,
            ),
            (
                SingleInterval(3, 5, Strand.PLUS),
                SingleInterval(0, 20, Strand.MINUS),
                DistanceType.OUTER,
                17,
            ),
            (
                SingleInterval(10, 20, Strand.MINUS),
                SingleInterval(10, 15, Strand.UNSTRANDED),
                DistanceType.OUTER,
                10,
            ),
            (
                SingleInterval(10, 20, Strand.MINUS),
                SingleInterval(13, 20, Strand.UNSTRANDED),
                DistanceType.OUTER,
                10,
            ),
            (
                SingleInterval(6, 8, Strand.PLUS),
                SingleInterval(0, 6, Strand.MINUS),
                DistanceType.STARTS,
                6,
            ),
            (
                SingleInterval(3, 5, Strand.PLUS),
                SingleInterval(0, 20, Strand.MINUS),
                DistanceType.STARTS,
                3,
            ),
            (
                SingleInterval(10, 20, Strand.MINUS),
                SingleInterval(10, 15, Strand.UNSTRANDED),
                DistanceType.STARTS,
                0,
            ),
            (
                SingleInterval(13, 20, Strand.UNSTRANDED),
                SingleInterval(0, 10, Strand.PLUS),
                DistanceType.ENDS,
                10,
            ),
            (
                SingleInterval(3, 5, Strand.PLUS),
                SingleInterval(0, 20, Strand.MINUS),
                DistanceType.ENDS,
                15,
            ),
            (
                SingleInterval(10, 20, Strand.MINUS),
                SingleInterval(13, 20, Strand.UNSTRANDED),
                DistanceType.ENDS,
                0,
            ),
            (
                SingleInterval(0, 5, Strand.PLUS),
                CompoundInterval.from_single_intervals([SingleInterval(7, 10, Strand.PLUS)]),
                DistanceType.INNER,
                2,
            ),
        ],
    )
    def test_distance_to(self, interval, other, distance_type, expected):
        print(interval, other, expected)  # %%%
        assert interval.distance_to(other, distance_type) == expected

    def test_distance_to_error(self):
        with pytest.raises(MismatchedParentException):
            SingleInterval(0, 1, Strand.PLUS, parent="seq1").distance_to(
                SingleInterval(0, 1, Strand.PLUS, parent="seq2")
            )

    @pytest.mark.parametrize(
        "interval,sequence_type,expected",
        [
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(id="parent", sequence_type="chr"),
                ),
                "chr",
                Parent(
                    id="parent",
                    sequence_type="chr",
                    strand=Strand.PLUS,
                    location=SingleInterval(5, 10, Strand.PLUS),
                ),
            ),
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="chr",
                        parent=Parent(id="grandparent", sequence_type="seq_chunk"),
                    ),
                ),
                "seq_chunk",
                Parent(id="grandparent", sequence_type="seq_chunk"),
            ),
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        parent=Parent(id="grandparent", sequence_type="seq_chunk"),
                    ),
                ),
                "seq_chunk",
                Parent(id="grandparent", sequence_type="seq_chunk"),
            ),
        ],
    )
    def test_first_ancestor_of_type(self, interval, sequence_type, expected):
        assert interval.first_ancestor_of_type(sequence_type) == expected

    @pytest.mark.parametrize(
        "interval,sequence_type",
        [
            # No parent
            (SingleInterval(5, 10, Strand.PLUS), "seq_chunk"),
            # Parent has no type
            (
                SingleInterval(5, 10, Strand.PLUS, parent="parent"),
                "seq_chunk",
            ),
            # Parent has wrong type, no additional ancestors
            (
                SingleInterval(5, 10, Strand.PLUS, parent=Parent(sequence_type="chr")),
                "seq_chunk",
            ),
            # Two ancestors with types, no matching type
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        sequence_type="chr",
                        parent="unknown",
                    ),
                ),
                "seq_chunk",
            ),
        ],
    )
    def test_first_ancestor_of_type_error(self, interval, sequence_type):
        with pytest.raises(NoSuchAncestorException):
            interval.first_ancestor_of_type(sequence_type)

    @pytest.mark.parametrize(
        "interval,sequence_type,expected",
        [
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(id="parent", sequence_type="chr"),
                ),
                "chr",
                True,
            ),
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="chr",
                        parent=Parent(id="grandparent", sequence_type="seq_chunk"),
                    ),
                ),
                "seq_chunk",
                True,
            ),
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        parent=Parent(id="grandparent", sequence_type="seq_chunk"),
                    ),
                ),
                "seq_chunk",
                True,
            ),
            # No parent
            (
                SingleInterval(5, 10, Strand.PLUS),
                "seq_chunk",
                False,
            ),
            # Parent has no type
            (
                SingleInterval(5, 10, Strand.PLUS, parent="parent"),
                "seq_chunk",
                False,
            ),
            # Parent has wrong type, no additional ancestors
            (
                SingleInterval(5, 10, Strand.PLUS, parent=Parent(sequence_type="chr")),
                "seq_chunk",
                False,
            ),
            # Two ancestors with types, no matching type
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        sequence_type="chr",
                        parent="unknown",
                    ),
                ),
                "seq_chunk",
                False,
            ),
        ],
    )
    def test_has_ancestor_of_type(self, interval, sequence_type, expected):
        assert interval.has_ancestor_of_type(sequence_type) is expected

    @pytest.mark.parametrize(
        "interval,sequence_type,expected",
        [
            # Parent
            (
                SingleInterval(
                    5,
                    10,
                    Strand.MINUS,
                    parent=Parent(sequence_type="chr"),
                ),
                "chr",
                SingleInterval(
                    5,
                    10,
                    Strand.MINUS,
                    parent=Parent(sequence_type="chr"),
                ),
            ),
            # Grandparent, parent has no type
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        parent=Parent(
                            id="grandparent",
                            strand=Strand.PLUS,
                            sequence_type="chr",
                            location=SingleInterval(100, 200, Strand.PLUS),
                        ),
                    ),
                ),
                "chr",
                SingleInterval(
                    105,
                    110,
                    Strand.PLUS,
                    Parent(
                        id="grandparent",
                        strand=Strand.PLUS,
                        sequence_type="chr",
                    ),
                ),
            ),
            # Grandparent, plus strand, parent has all attributes including its own parent
            (
                SingleInterval(
                    2,
                    4,
                    Strand.PLUS,
                    Parent(
                        id="parent",
                        sequence_type="chr",
                        strand=Strand.PLUS,
                        location=SingleInterval(2, 4, Strand.PLUS),
                        sequence=Sequence("AAAAA", Alphabet.NT_STRICT),
                        parent=Parent(
                            id="grandparent",
                            sequence_type="seq_chunk",
                            strand=Strand.MINUS,
                            location=SingleInterval(5, 10, Strand.MINUS),
                            parent="great_grandparent",
                        ),
                    ),
                ),
                "seq_chunk",
                SingleInterval(
                    6,
                    8,
                    Strand.MINUS,
                    Parent(
                        id="grandparent",
                        sequence_type="seq_chunk",
                        parent="great_grandparent",
                    ),
                ),
            ),
            # Grandparent, minus/plus
            (
                SingleInterval(
                    2,
                    4,
                    Strand.MINUS,
                    Parent(
                        id="parent",
                        sequence_type="seq_chunk",
                        parent=Parent(
                            id="grandparent",
                            sequence_type="unknown",
                            location=SingleInterval(3, 7, Strand.PLUS),
                        ),
                    ),
                ),
                "unknown",
                SingleInterval(
                    5,
                    7,
                    Strand.MINUS,
                    Parent(id="grandparent", sequence_type="unknown"),
                ),
            ),
            # Grandparent, minus/minus
            (
                SingleInterval(
                    2,
                    4,
                    Strand.MINUS,
                    Parent(
                        id="parent",
                        sequence_type="seq_chunk",
                        parent=Parent(
                            id="grandparent",
                            sequence_type="unknown",
                            strand=Strand.MINUS,
                            location=SingleInterval(3, 7, Strand.MINUS),
                        ),
                    ),
                ),
                "unknown",
                SingleInterval(
                    3,
                    5,
                    Strand.PLUS,
                    Parent(id="grandparent", sequence_type="unknown"),
                ),
            ),
            # Great grandparent, parent and grandparent have no type
            (
                SingleInterval(
                    2,
                    4,
                    Strand.PLUS,
                    Parent(
                        id="parent",
                        parent=Parent(
                            id="grandparent",
                            location=SingleInterval(0, 5, Strand.MINUS),
                            parent=Parent(
                                id="great_grandparent",
                                sequence_type="unknown",
                                location=SingleInterval(2, 8, Strand.PLUS),
                            ),
                        ),
                    ),
                ),
                "unknown",
                SingleInterval(
                    3,
                    5,
                    Strand.MINUS,
                    Parent(id="great_grandparent", sequence_type="unknown"),
                ),
            ),
            # Great grandparent, minus/minus
            (
                SingleInterval(
                    2,
                    4,
                    Strand.MINUS,
                    Parent(
                        id="parent",
                        sequence_type="seq_chunk",
                        parent=Parent(
                            id="grandparent",
                            location=SingleInterval(0, 5, Strand.PLUS),
                            parent=Parent(
                                id="great_grandparent",
                                sequence_type="unknown",
                                location=SingleInterval(2, 8, Strand.MINUS),
                            ),
                        ),
                    ),
                ),
                "unknown",
                SingleInterval(
                    4,
                    6,
                    Strand.PLUS,
                    Parent(id="great_grandparent", sequence_type="unknown"),
                ),
            ),
            # overlapping interval is child of a contiguous genome
            (
                CompoundInterval(
                    [5, 10],
                    [11, 20],
                    Strand.PLUS,
                    Parent(
                        id="chunk",
                        parent=Parent(
                            location=SingleInterval(
                                3,
                                100,
                                Strand.PLUS,
                                parent=Parent(id="genome", sequence_type="chromosome"),
                            )
                        ),
                    ),
                ),
                "chromosome",
                CompoundInterval(
                    [8, 13], [14, 23], Strand.PLUS, parent=Parent(id="genome", sequence_type="chromosome")
                ),
            ),
        ],
    )
    def test_lift_over_to_first_ancestor_of_type(self, interval, sequence_type, expected):
        assert interval.lift_over_to_first_ancestor_of_type(sequence_type) == expected

    @pytest.mark.parametrize(
        "interval,sequence_type,expected_error",
        [
            # No parent
            (
                SingleInterval(5, 10, Strand.PLUS),
                "seq_chunk",
                NoSuchAncestorException,
            ),
            # Parent has no type
            (
                SingleInterval(5, 10, Strand.PLUS, Parent(sequence_type="seq_chunk")),
                "unknown",
                NoSuchAncestorException,
            ),
            # Parent has wrong type, no additional ancestors
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    Parent(
                        sequence_type="seq_chunk",
                        parent="chr",
                    ),
                ),
                "unknown",
                NoSuchAncestorException,
            ),
            # Two ancestors with types, no matching type
            (
                SingleInterval(
                    2,
                    4,
                    Strand.MINUS,
                    Parent(
                        id="parent",
                        sequence_type="seq_chunk",
                        parent=Parent(id="grandparent", sequence_type="chr"),
                    ),
                ),
                "chr",
                NullParentException,
            ),
            # Parent has no location on grandparent
            (
                SingleInterval(
                    2,
                    4,
                    Strand.MINUS,
                    Parent(
                        id="parent",
                        sequence_type="seq_chunk",
                        parent=Parent(
                            id="grandparent",
                            sequence_type="unknown",
                            parent=Parent(
                                id="great_grandparent",
                                sequence_type="unknown",
                                location=SingleInterval(2, 8, Strand.MINUS),
                            ),
                        ),
                    ),
                ),
                "unknown",
                NullParentException,
            ),
        ],
    )
    def test_lift_over_to_first_ancestor_of_type_error(self, interval, sequence_type, expected_error):
        with pytest.raises(expected_error):
            interval.lift_over_to_first_ancestor_of_type(sequence_type)

    @pytest.mark.parametrize(
        "interval,sequence,expected",
        [
            # No parent
            (
                SingleInterval(0, 1, Strand.PLUS),
                Sequence("AA", Alphabet.NT_STRICT),
                False,
            ),
            # Has parent with given sequence
            (
                SingleInterval(
                    0,
                    1,
                    Strand.PLUS,
                    parent=Sequence("AA", Alphabet.NT_STRICT),
                ),
                Sequence("AA", Alphabet.NT_STRICT),
                True,
            ),
            # Has parent with sequence with one different attribute
            (
                SingleInterval(
                    0,
                    1,
                    Strand.PLUS,
                    parent=Sequence("AA", Alphabet.NT_STRICT),
                ),
                Sequence("AA", Alphabet.NT_STRICT, id="id"),
                False,
            ),
            # Grandparent
            (
                SingleInterval(
                    0,
                    1,
                    Strand.PLUS,
                    parent=Parent(
                        sequence=Sequence("AA", Alphabet.NT_STRICT),
                        parent=Sequence("AAT", Alphabet.NT_STRICT),
                    ),
                ),
                Sequence("AAT", Alphabet.NT_STRICT),
                True,
            ),
            # Has parent and grandparent with sequences, neither match
            (
                SingleInterval(
                    0,
                    1,
                    Strand.PLUS,
                    parent=Parent(
                        sequence=Sequence("AA", Alphabet.NT_STRICT),
                        parent=Sequence("AA", Alphabet.NT_STRICT),
                    ),
                ),
                Sequence("AAT", Alphabet.NT_STRICT),
                False,
            ),
        ],
    )
    def test_has_ancestor_sequence(self, interval, sequence, expected):
        assert interval.has_ancestor_sequence(sequence) == expected

    @pytest.mark.parametrize(
        "interval,sequence,expected",
        [
            # Parent
            (
                SingleInterval(
                    1,
                    2,
                    Strand.PLUS,
                    parent=Sequence("AAA", Alphabet.NT_STRICT),
                ),
                Sequence("AAA", Alphabet.NT_STRICT),
                SingleInterval(
                    1,
                    2,
                    Strand.PLUS,
                    parent=Sequence("AAA", Alphabet.NT_STRICT),
                ),
            ),
            # Grandparent, plus/minus
            (
                SingleInterval(
                    1,
                    2,
                    Strand.PLUS,
                    parent=Parent(
                        sequence=Sequence("AAA", Alphabet.NT_STRICT),
                        parent=Parent(
                            location=SingleInterval(2, 5, Strand.MINUS),
                            sequence=Sequence("TTTTT", Alphabet.NT_STRICT),
                        ),
                    ),
                ),
                Sequence("TTTTT", Alphabet.NT_STRICT),
                SingleInterval(
                    3,
                    4,
                    Strand.MINUS,
                    parent=Sequence("TTTTT", Alphabet.NT_STRICT),
                ),
            ),
        ],
    )
    def test_lift_over_to_sequence(self, interval, sequence, expected):
        assert interval.lift_over_to_sequence(sequence) == expected

    @pytest.mark.parametrize(
        "interval,sequence",
        [
            # No parent
            (SingleInterval(0, 1, Strand.PLUS), Sequence("AA", Alphabet.NT_STRICT)),
            # Has parent with sequence with different id
            (
                SingleInterval(
                    0,
                    1,
                    Strand.PLUS,
                    parent=Sequence("AA", Alphabet.NT_STRICT),
                ),
                Sequence("AA", Alphabet.NT_STRICT, id="id"),
            ),
            # Parent doesn't match, parent has empty parent
            (
                SingleInterval(
                    0,
                    1,
                    Strand.PLUS,
                    parent=Parent(sequence=Sequence("AA", Alphabet.NT_STRICT), parent=Parent()),
                ),
                Sequence("AA", Alphabet.NT_STRICT, id="id"),
            ),
        ],
    )
    def test_lift_over_to_sequence_error(self, interval, sequence):
        with pytest.raises(NoSuchAncestorException):
            interval.lift_over_to_sequence(sequence)

    @pytest.mark.parametrize(
        "interval1,interval2,expected",
        [
            (
                SingleInterval(
                    1,
                    1,
                    Strand.PLUS,
                    parent=Parent(id="parent", location=SingleInterval(1, 1, Strand.PLUS)),
                ),
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(id="parent", location=SingleInterval(5, 10, Strand.PLUS)),
                ),
                SingleInterval(5, 10, Strand.PLUS, parent="parent"),
            ),
            (
                SingleInterval(5, 10, Strand.MINUS),
                SingleInterval(20, 20, Strand.MINUS),
                SingleInterval(5, 10, Strand.MINUS),
            ),
            (
                SingleInterval(5, 10, Strand.UNSTRANDED),
                SingleInterval(15, 20, Strand.UNSTRANDED),
                CompoundInterval([5, 15], [10, 20], Strand.UNSTRANDED),
            ),
            (
                SingleInterval(15, 20, Strand.MINUS, parent="parent"),
                SingleInterval(5, 10, Strand.MINUS, parent="parent"),
                CompoundInterval([5, 15], [10, 20], Strand.MINUS, parent="parent"),
            ),
            (
                SingleInterval(5, 15, Strand.PLUS),
                SingleInterval(10, 20, Strand.PLUS),
                SingleInterval(5, 20, Strand.PLUS),
            ),
            (
                SingleInterval(10, 20, Strand.UNSTRANDED),
                SingleInterval(5, 15, Strand.UNSTRANDED),
                SingleInterval(5, 20, Strand.UNSTRANDED),
            ),
            (
                SingleInterval(0, 2, Strand.MINUS),
                SingleInterval(0, 2, Strand.MINUS),
                SingleInterval(0, 2, Strand.MINUS),
            ),
        ],
    )
    def test_union_single_interval(self, interval1, interval2, expected):
        assert interval1.union(interval2) == expected

    @pytest.mark.parametrize(
        "interval1,interval2,expected_exception",
        [
            (
                SingleInterval(5, 10, Strand.PLUS),
                SingleInterval(20, 30, Strand.MINUS),
                ValueError,
            ),
            (
                SingleInterval(5, 10, Strand.PLUS, parent="parent1"),
                SingleInterval(5, 10, Strand.PLUS, parent="parent2"),
                MismatchedParentException,
            ),
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(id="parent", sequence_type="chr"),
                ),
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(id="parent", sequence_type="unknown"),
                ),
                MismatchedParentException,
            ),
            (
                SingleInterval(5, 10, Strand.PLUS, parent="parent1"),
                SingleInterval(5, 10, Strand.PLUS),
                MismatchedParentException,
            ),
        ],
    )
    def test_union_single_interval_error(self, interval1, interval2, expected_exception):
        with pytest.raises(expected_exception):
            interval1.union(interval2)

    @pytest.mark.parametrize(
        "interval,compound_interval,expected",
        [
            # Empty single interval
            (
                SingleInterval(5, 5, Strand.PLUS),
                CompoundInterval([0, 10], [3, 13], Strand.PLUS),
                CompoundInterval([0, 10], [3, 13], Strand.PLUS),
            ),
            # Empty compound interval, both have parents
            (
                SingleInterval(
                    5,
                    10,
                    Strand.MINUS,
                    parent=Parent(id="parent", location=SingleInterval(5, 10, Strand.MINUS)),
                ),
                CompoundInterval(
                    [3],
                    [3],
                    Strand.MINUS,
                    parent=Parent(id="parent", location=CompoundInterval([3], [3], Strand.MINUS)),
                ),
                SingleInterval(5, 10, Strand.MINUS, parent="parent"),
            ),
            # Disjoint, one block, single interval on left
            (
                SingleInterval(5, 10, Strand.UNSTRANDED),
                CompoundInterval([20], [30], Strand.UNSTRANDED),
                CompoundInterval([5, 20], [10, 30], Strand.UNSTRANDED),
            ),
            # Disjoint, one block single interval on right
            (
                SingleInterval(20, 30, Strand.PLUS),
                CompoundInterval([5], [10], Strand.PLUS),
                CompoundInterval([5, 20], [10, 30], Strand.PLUS),
            ),
            # Disjoint, two blocks, both have parents, single interval on left
            (
                SingleInterval(
                    0,
                    3,
                    Strand.MINUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="unknown",
                        location=SingleInterval(0, 3, Strand.MINUS),
                    ),
                ),
                CompoundInterval(
                    [5, 15],
                    [10, 20],
                    Strand.MINUS,
                    parent=Parent(id="parent", sequence_type="unknown"),
                ),
                CompoundInterval(
                    [0, 5, 15],
                    [3, 10, 20],
                    Strand.MINUS,
                    parent=Parent(id="parent", sequence_type="unknown"),
                ),
            ),
            # Disjoint, two blocks, single interval on right
            (
                SingleInterval(20, 30, Strand.UNSTRANDED),
                CompoundInterval([0, 10], [5, 15], Strand.UNSTRANDED),
                CompoundInterval([0, 10, 20], [5, 15, 30], Strand.UNSTRANDED),
            ),
            # Single interval overlaps left block
            (
                SingleInterval(15, 25, Strand.PLUS),
                CompoundInterval([10, 30], [20, 40], Strand.PLUS),
                CompoundInterval([10, 30], [25, 40], Strand.PLUS),
            ),
            # Single interval overlaps middle block
            (
                SingleInterval(25, 35, Strand.MINUS),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.MINUS),
                CompoundInterval([10, 25, 50], [20, 40, 60], Strand.MINUS),
            ),
            # Single interval overlaps last block on the right
            (
                SingleInterval(55, 65, Strand.MINUS),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.MINUS),
                CompoundInterval([10, 30, 50], [20, 40, 65], Strand.MINUS),
            ),
            # Single interval contained in first block
            (
                SingleInterval(12, 13, Strand.PLUS),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.PLUS),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.PLUS),
            ),
            # Single interval equal to middle block
            (
                SingleInterval(30, 40, Strand.MINUS),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.MINUS),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.MINUS),
            ),
            # Single interval equal to an intron
            (
                SingleInterval(40, 50, Strand.UNSTRANDED),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.UNSTRANDED),
                CompoundInterval([10, 30], [20, 60], Strand.UNSTRANDED),
            ),
            # Single interval contained in an intron adjacent to a block
            (
                SingleInterval(20, 25, Strand.PLUS),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.PLUS),
                CompoundInterval([10, 30, 50], [25, 40, 60], Strand.PLUS),
            ),
            # Single interval contained in an intron not touching a block
            (
                SingleInterval(22, 25, Strand.MINUS),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.MINUS),
                CompoundInterval([10, 22, 30, 50], [20, 25, 40, 60], Strand.MINUS),
            ),
            # Single interval contains a block
            (
                SingleInterval(15, 35, Strand.UNSTRANDED),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.UNSTRANDED),
                CompoundInterval([10, 50], [40, 60], Strand.UNSTRANDED),
            ),
            # Single interval overlaps 3 blocks
            (
                SingleInterval(15, 60, Strand.PLUS),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.PLUS),
                SingleInterval(10, 60, Strand.PLUS),
            ),
            # Single interval overlaps two blocks
            (
                SingleInterval(35, 70, Strand.MINUS),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.MINUS),
                CompoundInterval([10, 30], [20, 70], Strand.MINUS),
            ),
            # Single interval contains entire compound interval
            (
                SingleInterval(0, 80, Strand.UNSTRANDED),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.UNSTRANDED),
                SingleInterval(0, 80, Strand.UNSTRANDED),
            ),
        ],
    )
    def test_union_compound_interval(self, interval, compound_interval, expected):
        assert interval.union(compound_interval) == expected

    @pytest.mark.parametrize(
        "interval,compound_interval,expected_exception",
        [
            # Opposite strands
            (
                SingleInterval(5, 10, Strand.PLUS),
                CompoundInterval([20], [30], Strand.MINUS),
                ValueError,
            ),
            # Different parents
            (
                SingleInterval(5, 10, Strand.PLUS, parent="parent1"),
                CompoundInterval([20], [30], Strand.PLUS, parent="parent2"),
                MismatchedParentException,
            ),
            # One has no parent
            (
                SingleInterval(5, 10, Strand.PLUS, parent="parent1"),
                CompoundInterval([20], [30], Strand.PLUS),
                MismatchedParentException,
            ),
        ],
    )
    def test_union_compound_interval_error(self, interval, compound_interval, expected_exception):
        with pytest.raises(expected_exception):
            interval.union(compound_interval)

    @pytest.mark.parametrize(
        "interval,other,exp",
        [
            # No overlap
            (
                SingleInterval(0, 5, Strand.PLUS),
                CompoundInterval([6, 10], [7, 11], Strand.PLUS),
                CompoundInterval([0, 6, 10], [5, 7, 11], Strand.PLUS),
            ),
            # Adjacent blocks
            (
                SingleInterval(0, 5, Strand.PLUS),
                CompoundInterval([5, 10], [6, 11], Strand.PLUS),
                CompoundInterval([0, 10], [6, 11], Strand.PLUS),
            ),
            # Overlapping blocks, has parent
            (
                SingleInterval(0, 5, Strand.PLUS, parent="seq"),
                SingleInterval(4, 10, Strand.PLUS, parent="seq"),
                CompoundInterval([0, 4], [5, 10], Strand.PLUS, parent="seq"),
            ),
        ],
    )
    def test_union_preserve_overlaps(self, interval, other, exp):
        assert interval.union_preserve_overlaps(other) == exp

    @pytest.mark.parametrize(
        "interval,other,exp_exception",
        [
            # Different parents
            (
                SingleInterval(0, 5, Strand.PLUS, parent="seq1"),
                SingleInterval(0, 5, Strand.PLUS, parent="seq2"),
                MismatchedParentException,
            ),
            # Different strands
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(0, 5, Strand.UNSTRANDED),
                InvalidStrandException,
            ),
        ],
    )
    def test_union_preserve_overlaps_error(self, interval, other, exp_exception):
        with pytest.raises(exp_exception):
            interval.union_preserve_overlaps(other)

    @pytest.mark.parametrize(
        "interval,other,match_strand,expected",
        [
            # Nested
            (
                SingleInterval(10, 20, Strand.MINUS),
                SingleInterval(15, 17, Strand.MINUS),
                True,
                SingleInterval(15, 17, Strand.MINUS),
            ),
            # Same span
            (
                SingleInterval(5, 10, Strand.UNSTRANDED),
                SingleInterval(5, 10, Strand.UNSTRANDED),
                True,
                SingleInterval(5, 10, Strand.UNSTRANDED),
            ),
            # Staggered; has parent
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(id="parent", location=SingleInterval(5, 10, Strand.PLUS)),
                ),
                SingleInterval(
                    9,
                    11,
                    Strand.PLUS,
                    parent=Parent(id="parent", location=SingleInterval(9, 11, Strand.PLUS)),
                ),
                True,
                SingleInterval(9, 10, Strand.PLUS, parent="parent"),
            ),
            # No overlap
            (
                SingleInterval(0, 1, Strand.PLUS),
                SingleInterval(1, 2, Strand.PLUS),
                True,
                EmptyLocation(),
            ),
            # Overlap but different strands, match_strand=False
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(0, 5, Strand.UNSTRANDED),
                False,
                SingleInterval(0, 5, Strand.PLUS),
            ),
            # Overlap but different strands
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(0, 5, Strand.UNSTRANDED),
                True,
                EmptyLocation(),
            ),
            # Overlap but different parents
            (
                SingleInterval(0, 5, Strand.PLUS, parent="parent1"),
                SingleInterval(0, 5, Strand.PLUS, parent="parent2"),
                True,
                EmptyLocation(),
            ),
        ],
    )
    def test_intersection_single_interval(self, interval, other, match_strand, expected):
        # full span has no effect when they are both single intervals
        assert (
            interval.intersection(other, match_strand)
            == interval.intersection(other, match_strand, full_span=True)
            == expected
        )
        if match_strand:
            assert (
                other.intersection(interval, match_strand)
                == interval.intersection(other, match_strand, full_span=True)
                == expected
            )

    def test_intersection_error(self):
        # Interval grandparents have different IDs
        seq1 = Sequence("A", Alphabet.NT_STRICT, parent=Sequence("A", Alphabet.NT_STRICT, id="ps1"))
        seq2 = Sequence("A", Alphabet.NT_STRICT, parent=Sequence("A", Alphabet.NT_STRICT, id="ps2"))
        with pytest.raises(MismatchedParentException):
            SingleInterval(0, 1, Strand.PLUS, parent=seq1).intersection(
                SingleInterval(0, 1, Strand.PLUS, parent=seq2), strict_parent_compare=True
            )

    @pytest.mark.parametrize(
        "compound_interval,single_interval,match_strand,expected",
        [
            # Different strands
            (
                CompoundInterval([5], [10], Strand.PLUS),
                SingleInterval(5, 10, Strand.MINUS),
                True,
                EmptyLocation(),
            ),
            # Different strands, match_strand=False
            (
                CompoundInterval([5], [11], Strand.PLUS),
                SingleInterval(5, 10, Strand.MINUS),
                False,
                SingleInterval(5, 10, Strand.MINUS),
            ),
            # No overlap
            (
                CompoundInterval([5], [10], Strand.PLUS),
                SingleInterval(10, 20, Strand.PLUS),
                True,
                EmptyLocation(),
            ),
            # Single interval contains compound interval
            (
                CompoundInterval([0, 5], [3, 10], Strand.MINUS),
                SingleInterval(0, 12, Strand.MINUS),
                True,
                CompoundInterval([0, 5], [3, 10], Strand.MINUS),
            ),
            # One block contains single interval
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                SingleInterval(3, 5, Strand.PLUS),
                True,
                SingleInterval(3, 5, Strand.PLUS),
            ),
            # Single interval overlaps multiple blocks
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                SingleInterval(3, 12, Strand.PLUS),
                True,
                CompoundInterval([3, 10], [5, 12], Strand.PLUS),
            ),
            # Both have parents and locations
            (
                CompoundInterval(
                    [5],
                    [10],
                    Strand.PLUS,
                    parent=Parent(id="parent", location=CompoundInterval([5], [10], Strand.PLUS)),
                ),
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(id="parent", location=SingleInterval(5, 10, Strand.PLUS)),
                ),
                True,
                SingleInterval(5, 10, Strand.PLUS, parent="parent"),
            ),
            # Different parents
            (
                CompoundInterval([5], [10], Strand.PLUS, parent="parent1"),
                SingleInterval(5, 10, Strand.PLUS, parent="parent2"),
                True,
                EmptyLocation(),
            ),
        ],
    )
    def test_intersection_compound_interval(self, compound_interval, single_interval, match_strand, expected):
        assert single_interval.intersection(compound_interval, match_strand=match_strand) == expected

    @pytest.mark.parametrize(
        "interval,other,match_strand,expected",
        [
            # Other completely contained in interval
            (
                SingleInterval(10, 20, Strand.MINUS),
                SingleInterval(12, 18, Strand.MINUS),
                True,
                CompoundInterval([10, 18], [12, 20], Strand.MINUS),
            ),
            # Different strand
            (
                SingleInterval(10, 20, Strand.MINUS),
                SingleInterval(12, 18, Strand.PLUS),
                True,
                SingleInterval(10, 20, Strand.MINUS),
            ),
            # Different strand, match_strand=False
            (
                SingleInterval(10, 20, Strand.MINUS),
                SingleInterval(12, 18, Strand.PLUS),
                False,
                CompoundInterval([10, 18], [12, 20], Strand.MINUS),
            ),
            # Same start, other contained in interval
            (
                SingleInterval(10, 20, Strand.PLUS),
                SingleInterval(10, 14, Strand.PLUS),
                True,
                SingleInterval(14, 20, Strand.PLUS),
            ),
            # Same end, other contained in interval
            (
                SingleInterval(10, 20, Strand.UNSTRANDED),
                SingleInterval(13, 20, Strand.UNSTRANDED),
                True,
                SingleInterval(10, 13, Strand.UNSTRANDED),
            ),
            # Overlap, interval left
            (
                SingleInterval(10, 20, Strand.PLUS),
                SingleInterval(15, 25, Strand.PLUS),
                True,
                SingleInterval(10, 15, Strand.PLUS),
            ),
            # Overlap, interval right, has parent
            (
                SingleInterval(
                    10,
                    20,
                    Strand.MINUS,
                    parent=Parent(id="parent", location=SingleInterval(10, 20, Strand.MINUS)),
                ),
                SingleInterval(
                    5,
                    15,
                    Strand.MINUS,
                    parent=Parent(id="parent", location=SingleInterval(5, 15, Strand.MINUS)),
                ),
                True,
                SingleInterval(15, 20, Strand.MINUS, parent="parent"),
            ),
            # Overlap, different parents
            (
                SingleInterval(10, 20, Strand.PLUS, parent="parent1"),
                SingleInterval(10, 12, Strand.PLUS, parent="parent2"),
                True,
                SingleInterval(10, 20, Strand.PLUS, parent="parent1"),
            ),
            # No overlap, interval left
            (
                SingleInterval(5, 7, Strand.PLUS),
                SingleInterval(7, 8, Strand.PLUS),
                True,
                SingleInterval(5, 7, Strand.PLUS),
            ),
            # No overlap, interval right
            (
                SingleInterval(5, 7, Strand.PLUS),
                SingleInterval(3, 5, Strand.PLUS),
                True,
                SingleInterval(5, 7, Strand.PLUS),
            ),
            # Interval completely contained in other, different parents
            (
                SingleInterval(10, 20, Strand.PLUS, parent="parent1"),
                SingleInterval(10, 25, Strand.PLUS, parent="parent2"),
                True,
                SingleInterval(10, 20, Strand.PLUS, parent="parent1"),
            ),
            # Interval completely contained in other, different strands
            (
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(0, 10, Strand.MINUS),
                True,
                SingleInterval(0, 5, Strand.PLUS),
            ),
            # Same span
            (
                SingleInterval(5, 10, Strand.PLUS),
                SingleInterval(5, 10, Strand.PLUS),
                True,
                EmptyLocation(),
            ),
            # Interval completely contained in other, same parent except location
            (
                SingleInterval(
                    5,
                    10,
                    Strand.PLUS,
                    parent=Parent(id="parent", location=SingleInterval(5, 10, Strand.PLUS)),
                ),
                SingleInterval(
                    4,
                    11,
                    Strand.PLUS,
                    parent=Parent(id="parent", location=SingleInterval(4, 11, Strand.PLUS)),
                ),
                True,
                EmptyLocation(),
            ),
        ],
    )
    def test_minus_single_interval(self, interval, other, match_strand, expected):
        assert interval.minus(other, match_strand) == expected

    @pytest.mark.parametrize(
        "location1,location2,match_strand,expected",
        [
            # No overlap
            (
                SingleInterval(0, 5, Strand.PLUS),
                CompoundInterval([5], [10], Strand.PLUS),
                True,
                SingleInterval(0, 5, Strand.PLUS),
            ),
            # Different strands
            (
                SingleInterval(0, 5, Strand.PLUS),
                CompoundInterval([0], [5], Strand.MINUS),
                True,
                SingleInterval(0, 5, Strand.PLUS),
            ),
            # Different parents
            (
                SingleInterval(0, 5, Strand.PLUS, parent="parent1"),
                CompoundInterval([0], [5], Strand.MINUS, parent="parent2"),
                True,
                SingleInterval(0, 5, Strand.PLUS, parent="parent1"),
            ),
        ],
    )
    def test_minus_compound_interval_no_overlap(self, location1, location2, match_strand, expected):
        assert location1.minus(location2, match_strand) == expected

    @pytest.mark.parametrize(
        "location1,location2,match_strand,expected",
        [
            # First contains second
            (
                SingleInterval(0, 20, Strand.MINUS),
                CompoundInterval([5, 10], [7, 12], Strand.MINUS),
                True,
                CompoundInterval([0, 7, 12], [5, 10, 20], Strand.MINUS),
            ),
            (
                SingleInterval(0, 20, Strand.MINUS),
                CompoundInterval([5, 10], [7, 20], Strand.MINUS),
                True,
                CompoundInterval([0, 7], [5, 10], Strand.MINUS),
            ),
            # Second contains first
            (
                SingleInterval(5, 10, Strand.PLUS),
                CompoundInterval([0, 20], [10, 30], Strand.PLUS),
                True,
                EmptyLocation(),
            ),
            # Multiple blocks overlap
            (
                SingleInterval(6, 11, Strand.MINUS),
                CompoundInterval([5, 10], [7, 12], Strand.MINUS),
                True,
                SingleInterval(7, 10, Strand.MINUS),
            ),
            (
                SingleInterval(6, 12, Strand.MINUS),
                CompoundInterval([5, 10], [7, 12], Strand.MINUS),
                True,
                SingleInterval(7, 10, Strand.MINUS),
            ),
            # Different strands, match_strand=False
            (
                SingleInterval(0, 20, Strand.PLUS),
                CompoundInterval([5, 10], [7, 12], Strand.MINUS),
                False,
                CompoundInterval([0, 7, 12], [5, 10, 20], Strand.PLUS),
            ),
        ],
    )
    def test_minus_compound_interval_with_overlap(self, location1, location2, match_strand, expected):
        assert location1.minus(location2, match_strand) == expected

    def test_minus_error(self):
        with pytest.raises(MismatchedParentException):
            SingleInterval(0, 1, Strand.PLUS, parent="seq1").minus(
                SingleInterval(0, 1, Strand.PLUS, parent="seq2"), strict_parent_compare=True
            )

    @pytest.mark.parametrize(
        "location,extend_left,extend_right,expected",
        [
            # Plus strand
            (SingleInterval(5, 10, Strand.PLUS), 5, 6, SingleInterval(0, 16, Strand.PLUS)),
            # Minus strand
            (SingleInterval(5, 10, Strand.MINUS), 5, 6, SingleInterval(0, 16, Strand.MINUS)),
            # Unstranded
            (SingleInterval(5, 10, Strand.UNSTRANDED), 5, 6, SingleInterval(0, 16, Strand.UNSTRANDED)),
            # Both distances zero
            (SingleInterval(5, 10, Strand.PLUS), 0, 0, SingleInterval(5, 10, Strand.PLUS)),
            # Has parent
            (
                SingleInterval(5, 10, Strand.MINUS, parent="parent"),
                5,
                6,
                SingleInterval(
                    0, 16, Strand.MINUS, parent=Parent(id="parent", location=SingleInterval(0, 16, Strand.MINUS))
                ),
            ),
        ],
    )
    def test_extend_absolute(self, location, extend_left, extend_right, expected):
        assert location.extend_absolute(extend_left, extend_right) == expected

    @pytest.mark.parametrize(
        "location,extend_left,extend_right,expected_exception",
        [
            # extend_left negative
            (SingleInterval(5, 10, Strand.PLUS), -1, 1, ValueError),
            # extend_right negative
            (SingleInterval(5, 10, Strand.PLUS), 1, -1, ValueError),
            # Extends left beyond parent
            (SingleInterval(5, 10, Strand.PLUS), 6, 1, InvalidPositionException),
            # Extends right beyond parent
            (
                SingleInterval(0, 5, Strand.PLUS, parent=Sequence("AAAAAA", Alphabet.NT_STRICT)),
                0,
                10,
                InvalidPositionException,
            ),
        ],
    )
    def test_extend_absolute_error(self, location, extend_left, extend_right, expected_exception):
        with pytest.raises(expected_exception):
            location.extend_absolute(extend_left, extend_right)

    @pytest.mark.parametrize(
        "location,extend_upstream,extend_downstream,expected",
        [
            # Plus strand
            (SingleInterval(5, 10, Strand.PLUS), 5, 6, SingleInterval(0, 16, Strand.PLUS)),
            # Minus strand
            (SingleInterval(5, 10, Strand.MINUS), 5, 3, SingleInterval(2, 15, Strand.MINUS)),
            # Both distances zero
            (SingleInterval(5, 10, Strand.PLUS), 0, 0, SingleInterval(5, 10, Strand.PLUS)),
            # Has parent
            (
                SingleInterval(5, 10, Strand.MINUS, parent="parent"),
                5,
                3,
                SingleInterval(
                    2, 15, Strand.MINUS, parent=Parent(id="parent", location=SingleInterval(2, 15, Strand.MINUS))
                ),
            ),
        ],
    )
    def test_extend_relative(self, location, extend_upstream, extend_downstream, expected):
        assert location.extend_relative(extend_upstream, extend_downstream) == expected

    @pytest.mark.parametrize(
        "location,extend_upstream,extend_downstream,expected_exception",
        [
            # Unstranded
            (SingleInterval(5, 10, Strand.UNSTRANDED), 1, 1, InvalidStrandException),
            # extend_left negative
            (SingleInterval(5, 10, Strand.PLUS), -1, 1, ValueError),
            # extend_right negative
            (SingleInterval(5, 10, Strand.PLUS), 1, -1, ValueError),
            # Extends left beyond parent
            (SingleInterval(5, 10, Strand.MINUS), 0, 6, InvalidPositionException),
            (SingleInterval(5, 10, Strand.PLUS), 6, 0, InvalidPositionException),
            # Extends right beyond parent
            (
                SingleInterval(0, 5, Strand.PLUS, parent=Sequence("AAAAAA", Alphabet.NT_STRICT)),
                0,
                10,
                InvalidPositionException,
            ),
            (
                SingleInterval(0, 5, Strand.MINUS, parent=Sequence("AAAAAA", Alphabet.NT_STRICT)),
                10,
                0,
                InvalidPositionException,
            ),
        ],
    )
    def test_extend_relative_error(self, location, extend_upstream, extend_downstream, expected_exception):
        with pytest.raises(expected_exception):
            location.extend_relative(extend_upstream, extend_downstream)

    @pytest.mark.parametrize(
        "interval,other,match_strand,expected",
        [
            # Contains, same strand, toggle match_strand
            (
                SingleInterval(0, 10, Strand.PLUS),
                SingleInterval(5, 7, Strand.PLUS),
                True,
                True,
            ),
            (
                SingleInterval(0, 10, Strand.PLUS),
                SingleInterval(5, 7, Strand.PLUS),
                False,
                True,
            ),
            # Contains, different strand, toggle match_strand
            (
                SingleInterval(0, 10, Strand.PLUS),
                SingleInterval(5, 7, Strand.MINUS),
                True,
                False,
            ),
            (
                SingleInterval(0, 10, Strand.PLUS),
                SingleInterval(5, 7, Strand.MINUS),
                False,
                True,
            ),
            # Equal, same strand
            (
                SingleInterval(5, 10, Strand.UNSTRANDED),
                SingleInterval(5, 10, Strand.UNSTRANDED),
                True,
                True,
            ),
            # Overlapping
            (
                SingleInterval(5, 10, Strand.PLUS),
                SingleInterval(6, 11, Strand.PLUS),
                True,
                False,
            ),
            # No overlap
            (
                SingleInterval(5, 6, Strand.PLUS),
                SingleInterval(6, 10, Strand.PLUS),
                True,
                False,
            ),
            # Contains, same location, different parents
            (
                SingleInterval(5, 10, Strand.PLUS, parent="parent1"),
                SingleInterval(5, 10, Strand.PLUS, parent="parent2"),
                True,
                False,
            ),
            # Contains, has parent with location
            (
                SingleInterval(
                    0,
                    10,
                    Strand.PLUS,
                    parent=Parent(id="parent", location=SingleInterval(0, 10, Strand.PLUS)),
                ),
                SingleInterval(
                    5,
                    7,
                    Strand.PLUS,
                    parent=Parent(id="parent", location=SingleInterval(5, 7, Strand.PLUS)),
                ),
                True,
                True,
            ),
            # Contained in other
            (
                SingleInterval(5, 7, Strand.PLUS),
                SingleInterval(5, 8, Strand.PLUS),
                True,
                False,
            ),
        ],
    )
    def test_contains_single_interval(self, interval, other, match_strand, expected):
        assert (
            interval.contains(other, match_strand) == interval.contains(other, match_strand, full_span=True) == expected
        )

    @pytest.mark.parametrize(
        "interval,other,match_strand,expected",
        [
            # Contains, different strand, toggle match_strand
            (
                SingleInterval(0, 20, Strand.PLUS),
                CompoundInterval([3, 10], [5, 12], Strand.MINUS),
                True,
                False,
            ),
            (
                SingleInterval(0, 20, Strand.PLUS),
                CompoundInterval([3, 10], [5, 12], Strand.MINUS),
                False,
                True,
            ),
            # Overlapping
            (
                SingleInterval(0, 10, Strand.PLUS),
                CompoundInterval([5, 20], [15, 25], Strand.PLUS),
                False,
                False,
            ),
            # No overlap
            (
                SingleInterval(30, 40, Strand.PLUS),
                CompoundInterval([5, 20], [15, 25], Strand.PLUS),
                False,
                False,
            ),
            # Contains one block of other
            (
                SingleInterval(0, 20, Strand.PLUS),
                CompoundInterval([5, 20], [15, 25], Strand.PLUS),
                False,
                False,
            ),
            # Contained in other
            (
                SingleInterval(5, 10, Strand.PLUS),
                CompoundInterval([5, 20], [15, 25], Strand.PLUS),
                False,
                False,
            ),
            # different parent
            (
                SingleInterval(0, 20, Strand.PLUS, parent="parent1"),
                CompoundInterval([3, 10], [5, 12], Strand.MINUS, parent="parent2"),
                False,
                False,
            ),
        ],
    )
    def test_contains_compound_interval(self, interval, other, match_strand, expected):
        assert interval.contains(other, match_strand) is expected

    @pytest.mark.parametrize(
        "interval,other,match_strand,expected",
        [
            # Overlapping
            (
                SingleInterval(0, 30, Strand.PLUS),
                CompoundInterval([5, 20], [15, 25], Strand.PLUS),
                False,
                True,
            ),
            # Overlapping but different parent
            (
                SingleInterval(0, 30, Strand.PLUS, parent="parent1"),
                CompoundInterval([5, 20], [15, 25], Strand.PLUS, parent="parent2"),
                False,
                False,
            ),
            # No overlap
            (
                SingleInterval(30, 40, Strand.PLUS),
                CompoundInterval([5, 20], [15, 25], Strand.PLUS),
                False,
                False,
            ),
            # Contains one block of other
            (
                SingleInterval(0, 20, Strand.PLUS),
                CompoundInterval([5, 20], [15, 25], Strand.PLUS),
                False,
                False,
            ),
            # Contained in other
            (
                SingleInterval(5, 10, Strand.PLUS),
                CompoundInterval([5, 20], [15, 25], Strand.PLUS),
                False,
                False,
            ),
            # Overlapping
            (
                SingleInterval(0, 30, Strand.PLUS),
                CompoundInterval([5, 20], [15, 25], Strand.MINUS),
                False,
                True,
            ),
            # No overlap
            (
                SingleInterval(30, 40, Strand.PLUS),
                CompoundInterval([5, 20], [15, 25], Strand.MINUS),
                False,
                False,
            ),
            # Contains one block of other
            (
                SingleInterval(0, 20, Strand.PLUS),
                CompoundInterval([5, 20], [15, 25], Strand.MINUS),
                False,
                False,
            ),
            # Contained in other
            (
                SingleInterval(5, 10, Strand.PLUS),
                CompoundInterval([5, 20], [15, 25], Strand.MINUS),
                False,
                False,
            ),
        ],
    )
    def test_contains_compound_interval_full_span(self, interval, other, match_strand, expected):
        assert interval.contains(other, match_strand, full_span=True) is expected

    def test_contains_error(self):
        with pytest.raises(MismatchedParentException):
            SingleInterval(0, 1, Strand.PLUS, parent="seq1").contains(
                SingleInterval(0, 1, Strand.PLUS, parent="seq2"), strict_parent_compare=True
            )

    @pytest.mark.parametrize(
        "location,expected_output",
        [
            (SingleInterval(0, 10, Strand.PLUS), FeatureLocation(ExactPosition(0), ExactPosition(10), strand=1)),
            (SingleInterval(0, 10, Strand.MINUS), FeatureLocation(ExactPosition(0), ExactPosition(10), strand=-1)),
        ],
    )
    def test_biopython(self, location, expected_output):
        assert location.to_feature_location() == expected_output
