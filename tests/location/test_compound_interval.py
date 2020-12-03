import pytest

from Bio.SeqFeature import FeatureLocation, CompoundLocation, ExactPosition

from inscripta.biocantor.exc import (
    NoSuchAncestorException,
    InvalidStrandException,
    InvalidPositionException,
)
from inscripta.biocantor.location.distance import DistanceType
from inscripta.biocantor.location.location_impl import (
    SingleInterval,
    CompoundInterval,
    EmptyLocation,
)
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.parent import Parent
from inscripta.biocantor.sequence import Sequence


class TestCompoundInterval:
    @pytest.mark.parametrize(
        "starts,ends,parent,expected_exception",
        [
            # Lists of starts and ends are empty
            ([], [], None, ValueError),
            # Different numbers of start and end positions
            ([0], [3, 5], None, ValueError),
            # An end is less than a start
            ([10], [9], None, InvalidPositionException),
            # Total length of intervals is longer than parent
            ([0], [3], Parent(sequence=Sequence("AA", Alphabet.NT_STRICT)), ValueError),
        ],
    )
    def test_init_error(self, starts, ends, parent, expected_exception):
        with pytest.raises(expected_exception):
            CompoundInterval(starts, ends, Strand.PLUS, parent)

    @pytest.mark.parametrize(
        "intervals,expected_exception",
        [
            # List of intervals is empty
            ([], ValueError),
            # Intervals have different parents
            (
                [
                    SingleInterval(0, 5, Strand.PLUS, parent="parent1"),
                    SingleInterval(6, 7, Strand.PLUS, parent="parent2"),
                ],
                ValueError,
            ),
            # Different strands
            (
                [SingleInterval(0, 5, Strand.PLUS), SingleInterval(6, 7, Strand.MINUS)],
                ValueError,
            ),
        ],
    )
    def test_from_single_intervals_error(self, intervals, expected_exception):
        with pytest.raises(expected_exception):
            CompoundInterval.from_single_intervals(intervals)

    @pytest.mark.parametrize(
        "location,expected",
        [
            # One block
            (CompoundInterval([3], [5], Strand.PLUS), 1),
            # Two adjacent blocks
            (CompoundInterval([3, 5], [5, 7], Strand.PLUS), 2),
            # Three blocks
            (CompoundInterval([3, 5, 10], [5, 7, 15], Strand.PLUS), 3),
        ],
    )
    def test_len(self, location, expected):
        assert location.num_blocks == expected

    def test_parent(self):
        # Parent gets location from self
        parent = Parent(id="parent")
        location = CompoundInterval([1, 5], [3, 7], Strand.PLUS, parent=parent)
        assert location.parent == Parent(id="parent", location=CompoundInterval([1, 5], [3, 7], Strand.PLUS))

    def test_parent_id(self):
        parent = Parent(id="parent")
        location = CompoundInterval([1, 5], [3, 7], Strand.PLUS, parent=parent)
        assert location.parent_id == "parent"

    def test_parent_type(self):
        parent = Parent(id="parent", sequence_type="chr")
        location = CompoundInterval([1, 5], [3, 7], Strand.PLUS, parent=parent)
        assert location.parent_type == "chr"

    def test_strand(self):
        assert CompoundInterval([5], [6], Strand.UNSTRANDED).strand == Strand.UNSTRANDED

    def test_start(self):
        assert CompoundInterval([5, 10], [6, 15], Strand.UNSTRANDED).start == 5

    def test_end(self):
        assert CompoundInterval([5, 10], [6, 15], Strand.UNSTRANDED).end == 15

    def test_blocks(self):
        # Regular constructor: intervals get sorted
        assert CompoundInterval([10, 5], [15, 6], Strand.PLUS).blocks == [
            SingleInterval(5, 6, Strand.PLUS),
            SingleInterval(10, 15, Strand.PLUS),
        ]
        # From single intervals: intervals get sorted
        assert CompoundInterval.from_single_intervals(
            [SingleInterval(10, 15, Strand.PLUS), SingleInterval(5, 6, Strand.PLUS)]
        ).blocks == [
            SingleInterval(5, 6, Strand.PLUS),
            SingleInterval(10, 15, Strand.PLUS),
        ]

    def test_overlapping(self):
        overlapping = CompoundInterval([0, 3], [5, 7], Strand.PLUS)
        assert overlapping.blocks == [
            SingleInterval(0, 5, Strand.PLUS),
            SingleInterval(3, 7, Strand.PLUS),
        ]
        assert overlapping.is_overlapping
        # not overlapping
        normal = CompoundInterval([10, 5], [15, 6], Strand.PLUS)
        assert not normal.is_overlapping

    def test_scan_blocks(self):
        assert list(CompoundInterval([10, 5, 20], [15, 6, 25], Strand.PLUS).scan_blocks()) == [
            SingleInterval(5, 6, Strand.PLUS),
            SingleInterval(10, 15, Strand.PLUS),
            SingleInterval(20, 25, Strand.PLUS),
        ]
        assert list(CompoundInterval([10, 5, 20], [15, 6, 25], Strand.MINUS).scan_blocks()) == [
            SingleInterval(20, 25, Strand.MINUS),
            SingleInterval(10, 15, Strand.MINUS),
            SingleInterval(5, 6, Strand.MINUS),
        ]

    def test_scan_blocks_error(self):
        with pytest.raises(InvalidStrandException):
            list(CompoundInterval([10, 5], [15, 6], Strand.UNSTRANDED).scan_blocks())

    @pytest.mark.parametrize(
        "location,expected",
        [
            # One interval
            (CompoundInterval([5], [10], Strand.PLUS), True),
            # Two adjacent intervals
            (CompoundInterval([5, 10], [10, 15], Strand.PLUS), True),
            # Two non-adjacent intervals
            (CompoundInterval([5, 10], [6, 11], Strand.PLUS), False),
        ],
    )
    def test_is_contiguous(self, location, expected):
        assert location.is_contiguous is expected

    @pytest.mark.parametrize(
        "object1,object2,expected",
        [
            # A SingleInterval and a CompoundInterval made from the same one block
            (
                SingleInterval(0, 1, Strand.PLUS),
                CompoundInterval([0], [1], Strand.PLUS),
                False,
            ),
            # One interval has different start
            (
                CompoundInterval([0, 5], [3, 7], Strand.PLUS),
                CompoundInterval([0, 4], [3, 7], Strand.PLUS),
                False,
            ),
            # Different strand
            (
                CompoundInterval([0, 5], [3, 7], Strand.PLUS),
                CompoundInterval([0, 5], [3, 7], Strand.MINUS),
                False,
            ),
            # Different parent
            (
                CompoundInterval([0, 5], [3, 7], Strand.PLUS, parent="parent1"),
                CompoundInterval([0, 5], [3, 7], Strand.PLUS, parent="parent2"),
                False,
            ),
            # Equal
            (
                CompoundInterval([0, 5], [3, 7], Strand.PLUS, parent="parent"),
                CompoundInterval([0, 5], [3, 7], Strand.PLUS, parent="parent"),
                True,
            ),
        ],
    )
    def test_equals(self, object1, object2, expected):
        assert (object1 == object2) is expected

    @pytest.mark.parametrize(
        "location,expected",
        [
            # One interval minus strand
            (
                CompoundInterval(
                    [0],
                    [3],
                    Strand.MINUS,
                    parent=Parent(sequence=Sequence("CCAA", Alphabet.NT_STRICT)),
                ),
                Sequence("TGG", Alphabet.NT_STRICT),
            ),
            # Two adjacent intervals plus strand
            (
                CompoundInterval(
                    [0, 2],
                    [2, 4],
                    Strand.PLUS,
                    parent=Sequence("AACCGGTT", Alphabet.NT_STRICT),
                ),
                Sequence("AACC", Alphabet.NT_STRICT),
            ),
            # Two non-adjacent intervals; minus strand; parent has all info; parent gets stripped
            (
                CompoundInterval(
                    [0, 3],
                    [2, 5],
                    Strand.MINUS,
                    parent=Parent(
                        sequence=Sequence("AACCGGTT", Alphabet.NT_STRICT),
                        location=CompoundInterval([0, 3], [2, 5], Strand.MINUS),
                        sequence_type="chr",
                    ),
                ),
                Sequence("CGTT", Alphabet.NT_STRICT),
            ),
        ],
    )
    def test_extract_sequence(self, location, expected):
        assert location.extract_sequence() == expected

    @pytest.mark.parametrize(
        "location,expected_exception",
        [
            # Unstranded
            (
                CompoundInterval(
                    [0],
                    [1],
                    Strand.UNSTRANDED,
                    parent=Sequence("AAA", Alphabet.NT_STRICT),
                ),
                InvalidStrandException,
            ),
            # No parent
            (CompoundInterval([0], [1], Strand.PLUS), ValueError),
            # Parent has no sequence
            (
                CompoundInterval([0], [1], Strand.PLUS, parent=Parent(id="parent")),
                ValueError,
            ),
        ],
    )
    def test_extract_sequence_error(self, location, expected_exception):
        with pytest.raises(expected_exception):
            location.extract_sequence()

    @pytest.mark.parametrize(
        "location,parent_pos,expected",
        [
            # In first intron; minus strand
            (CompoundInterval([0, 6, 20], [3, 8, 25], Strand.MINUS), 2, 7),
            # In first intron; plus strand
            (CompoundInterval([0, 6, 20], [3, 8, 25], Strand.PLUS), 2, 2),
            # In last intron; plus strand
            (CompoundInterval([0, 6, 20], [3, 8, 25], Strand.PLUS), 24, 9),
            # Last position before an intron; minus strand
            (CompoundInterval([0, 6, 20], [3, 8, 25], Strand.MINUS), 7, 5),
        ],
    )
    def test_parent_to_relative_pos(self, location, parent_pos, expected):
        assert location.parent_to_relative_pos(parent_pos) == expected

    @pytest.mark.parametrize(
        "location,parent_pos,expected_exception",
        [
            # Position is negative
            (
                CompoundInterval([0, 6, 20], [3, 8, 25], Strand.MINUS),
                -1,
                InvalidPositionException,
            ),
            # Position is greater than length of parent
            (
                CompoundInterval(
                    [0],
                    [1],
                    Strand.PLUS,
                    parent=Sequence("AAA", Alphabet.NT_STRICT),
                ),
                5,
                InvalidPositionException,
            ),
            # Position is adjacent to location
            (
                CompoundInterval([0, 6, 20], [3, 8, 25], Strand.MINUS),
                25,
                InvalidPositionException,
            ),
            # Position is in "intron"
            (
                CompoundInterval([0, 6, 20], [3, 8, 25], Strand.MINUS),
                10,
                InvalidPositionException,
            ),
            # Unstranded
            (
                CompoundInterval([0, 6, 20], [3, 8, 25], Strand.UNSTRANDED),
                3,
                InvalidStrandException,
            ),
        ],
    )
    def test_parent_to_relative_pos_error(self, location, parent_pos, expected_exception):
        with pytest.raises(expected_exception):
            location.parent_to_relative_pos(parent_pos)

    @pytest.mark.parametrize(
        "location,relative_pos,expected",
        [
            # In first intron; minus strand
            (CompoundInterval([0, 6, 20], [3, 8, 25], Strand.MINUS), 7, 2),
            # In first intron; plus strand
            (CompoundInterval([0, 6, 20], [3, 8, 25], Strand.PLUS), 2, 2),
            # In last intron; plus strand
            (CompoundInterval([0, 6, 20], [3, 8, 25], Strand.PLUS), 9, 24),
            # Last position before an intron; minus strand
            (CompoundInterval([0, 6, 20], [3, 8, 25], Strand.MINUS), 5, 7),
        ],
    )
    def test_relative_to_parent_pos(self, location, relative_pos, expected):
        assert location.relative_to_parent_pos(relative_pos) == expected

    @pytest.mark.parametrize(
        "location,relative_pos,expected_exception",
        [
            # Position is negative
            (CompoundInterval([0], [2], Strand.PLUS), -1, InvalidPositionException),
            # Position is greater than length of location
            (CompoundInterval([0], [2], Strand.PLUS), 2, InvalidPositionException),
            # Unstranded
            (CompoundInterval([0], [2], Strand.UNSTRANDED), 1, InvalidStrandException),
        ],
    )
    def test_relative_to_parent_pos_error(self, location, relative_pos, expected_exception):
        with pytest.raises(expected_exception):
            location.relative_to_parent_pos(relative_pos)

    def test_parent_to_relative_location_empty_location(self):
        assert (
            CompoundInterval([0], [10], Strand.PLUS, parent="parent").parent_to_relative_location(EmptyLocation())
            is EmptyLocation()
        )

    @pytest.mark.parametrize(
        "location,parent_location,expected",
        [
            # Both have parent
            (
                CompoundInterval([0], [10], Strand.PLUS, parent="parent"),
                SingleInterval(5, 7, Strand.MINUS, parent="parent"),
                SingleInterval(5, 7, Strand.MINUS),
            ),
            # Parent location is unstranded
            (
                CompoundInterval([0, 10], [5, 15], Strand.MINUS),
                SingleInterval(0, 12, Strand.UNSTRANDED),
                SingleInterval(3, 10, Strand.UNSTRANDED),
            ),
            # Parent location contains entire location; same strand
            (
                CompoundInterval([1, 10], [2, 14], Strand.MINUS),
                SingleInterval(0, 25, Strand.MINUS),
                SingleInterval(0, 5, Strand.PLUS),
            ),
            # Parent location overlaps two blocks; opposite strands
            (
                CompoundInterval([0, 10, 20], [3, 13, 23], Strand.PLUS),
                SingleInterval(7, 22, Strand.MINUS),
                SingleInterval(3, 8, Strand.MINUS),
            ),
            # Parent location equals location
            (
                CompoundInterval([3], [5], Strand.PLUS),
                SingleInterval(3, 5, Strand.PLUS),
                SingleInterval(0, 2, Strand.PLUS),
            ),
        ],
    )
    def test_parent_to_relative_location_single_interval(self, location, parent_location, expected):
        assert location.parent_to_relative_location(parent_location) == expected

    @pytest.mark.parametrize(
        "location,parent_location,expected",
        [
            # Both have parent
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS, parent="parent"),
                CompoundInterval([12, 20], [18, 25], Strand.MINUS, parent="parent"),
                SingleInterval(7, 10, Strand.MINUS),
            ),
            # Parent location is unstranded; overlaps two blocks
            (
                CompoundInterval([0, 10], [7, 17], Strand.PLUS),
                CompoundInterval([5, 20], [12, 30], Strand.UNSTRANDED),
                SingleInterval(5, 9, Strand.UNSTRANDED),
            ),
            # Parent location contains entire location; same strand
            (
                CompoundInterval([5, 10, 15], [7, 13, 19], Strand.MINUS),
                CompoundInterval([3, 8, 14], [8, 13, 20], Strand.MINUS),
                SingleInterval(0, 9, Strand.PLUS),
            ),
            # Parent location equals location
            (
                CompoundInterval([0, 3], [1, 10], Strand.PLUS),
                CompoundInterval([0, 3], [1, 10], Strand.PLUS),
                SingleInterval(0, 8, Strand.PLUS),
            ),
            # Parent location has multiple blocks within multiple blocks of location
            (
                CompoundInterval([0, 20, 40], [10, 30, 50], Strand.MINUS),
                CompoundInterval([3, 24, 41], [22, 27, 45], Strand.PLUS),
                CompoundInterval([5, 13, 18], [9, 16, 27], Strand.MINUS),
            ),
        ],
    )
    def test_parent_to_relative_location_compound_interval(self, location, parent_location, expected):
        assert location.parent_to_relative_location(parent_location) == expected

    @pytest.mark.parametrize(
        "location,parent_location,expected_exception",
        [
            # Different parents
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS, parent="parent1"),
                CompoundInterval([12, 20], [18, 25], Strand.MINUS, parent="parent2"),
                ValueError,
            ),
            # Location on parent is unstranded
            (
                CompoundInterval([3], [5], Strand.UNSTRANDED),
                SingleInterval(3, 5, Strand.PLUS),
                InvalidStrandException,
            ),
            # Parent location is outside location
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                SingleInterval(20, 25, Strand.PLUS),
                ValueError,
            ),
            # Parent location is in "intron"
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                SingleInterval(7, 8, Strand.PLUS),
                ValueError,
            ),
        ],
    )
    def test_parent_to_relative_location_error(self, location, parent_location, expected_exception):
        with pytest.raises(expected_exception):
            location.parent_to_relative_location(parent_location)

    @pytest.mark.parametrize(
        "location,start,end,strand,expected",
        [
            # Unstranded; location has parent
            (
                CompoundInterval(
                    [10, 30],
                    [13, 40],
                    Strand.MINUS,
                    parent=Parent(
                        id="parent",
                        location=CompoundInterval([10, 30], [13, 40], Strand.MINUS),
                    ),
                ),
                7,
                12,
                Strand.UNSTRANDED,
                CompoundInterval([11, 30], [13, 33], Strand.UNSTRANDED, parent="parent"),
            ),
            # Entire length of location; opposite strands
            (
                CompoundInterval([10, 15], [12, 17], Strand.PLUS),
                0,
                4,
                Strand.MINUS,
                CompoundInterval([10, 15], [12, 17], Strand.MINUS),
            ),
            # Spanning multiple blocks; both minus strand
            (
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.MINUS),
                5,
                25,
                Strand.MINUS,
                CompoundInterval([15, 30, 50], [20, 40, 55], Strand.PLUS),
            ),
            # Start == end
            (CompoundInterval([10], [20], Strand.PLUS), 5, 5, Strand.PLUS, SingleInterval(15, 15, Strand.PLUS)),
            # Overlapping blocks
            (
                CompoundInterval([0, 9], [10, 20], Strand.PLUS),
                8,
                12,
                Strand.PLUS,
                CompoundInterval([8, 9], [10, 11], Strand.PLUS),
            ),
        ],
    )
    def test_relative_interval_to_parent_location(self, location, start, end, strand, expected):
        assert location.relative_interval_to_parent_location(start, end, strand) == expected

    @pytest.mark.parametrize(
        "location,start,end,strand,expected_exception",
        [
            # Location on parent is unstranded
            (
                CompoundInterval([10, 15], [12, 17], Strand.UNSTRANDED),
                0,
                4,
                Strand.MINUS,
                InvalidStrandException,
            ),
            # Start is negative
            (
                CompoundInterval([10, 15], [12, 17], Strand.PLUS),
                -1,
                4,
                Strand.MINUS,
                InvalidPositionException,
            ),
            # End is greater than length of location
            (
                CompoundInterval([10, 15], [12, 17], Strand.PLUS),
                0,
                5,
                Strand.MINUS,
                InvalidPositionException,
            ),
            # Start > end
            (
                CompoundInterval([10, 15], [12, 17], Strand.PLUS),
                2,
                1,
                Strand.MINUS,
                InvalidPositionException,
            ),
        ],
    )
    def test_relative_interval_to_parent_location_error(self, location, start, end, strand, expected_exception):
        with pytest.raises(expected_exception):
            location.relative_interval_to_parent_location(start, end, strand)

    @pytest.mark.parametrize(
        "location,window_size,step_size,start_pos,expected",
        [
            # Step size greater than window size, windows go to end of location
            (
                CompoundInterval([0, 20, 40], [10, 27, 53], Strand.PLUS),
                2,
                5,
                3,
                [
                    SingleInterval(3, 5, Strand.PLUS),
                    SingleInterval(8, 10, Strand.PLUS),
                    SingleInterval(23, 25, Strand.PLUS),
                    SingleInterval(41, 43, Strand.PLUS),
                    SingleInterval(46, 48, Strand.PLUS),
                    SingleInterval(51, 53, Strand.PLUS),
                ],
            ),
            # Start pos is in second exon; minus strand
            (
                CompoundInterval([0, 20, 40], [10, 27, 54], Strand.MINUS),
                5,
                3,
                18,
                [
                    CompoundInterval([8, 20], [10, 23], Strand.MINUS),
                    SingleInterval(5, 10, Strand.MINUS),
                    SingleInterval(2, 7, Strand.MINUS),
                ],
            ),
            # Step size = 1; plus strand
            (
                CompoundInterval([0, 7], [3, 10], Strand.PLUS),
                2,
                1,
                0,
                [
                    SingleInterval(0, 2, Strand.PLUS),
                    SingleInterval(1, 3, Strand.PLUS),
                    CompoundInterval([2, 7], [3, 8], Strand.PLUS),
                    SingleInterval(7, 9, Strand.PLUS),
                    SingleInterval(8, 10, Strand.PLUS),
                ],
            ),
            # Step size = 2; minus strand
            (
                CompoundInterval([10, 30], [20, 37], Strand.MINUS),
                5,
                2,
                4,
                [
                    CompoundInterval([18, 30], [20, 33], Strand.MINUS),
                    CompoundInterval([16, 30], [20, 31], Strand.MINUS),
                    SingleInterval(14, 19, Strand.MINUS),
                    SingleInterval(12, 17, Strand.MINUS),
                    SingleInterval(10, 15, Strand.MINUS),
                ],
            ),
            # Window size = 1, step size = 1; minus strand
            (
                CompoundInterval([3, 10], [5, 13], Strand.MINUS),
                1,
                1,
                0,
                [
                    SingleInterval(12, 13, Strand.MINUS),
                    SingleInterval(11, 12, Strand.MINUS),
                    SingleInterval(10, 11, Strand.MINUS),
                    SingleInterval(4, 5, Strand.MINUS),
                    SingleInterval(3, 4, Strand.MINUS),
                ],
            ),
            # Window size is equal to location length
            (
                CompoundInterval([3, 10], [5, 13], Strand.PLUS),
                5,
                1,
                0,
                [CompoundInterval([3, 10], [5, 13], Strand.PLUS)],
            ),
            # Start pos is last position of location, window size 1
            (
                CompoundInterval([3, 10], [5, 13], Strand.PLUS),
                1,
                1,
                4,
                [SingleInterval(12, 13, Strand.PLUS)],
            ),
            # Overlapping blocks
            (
                CompoundInterval([0, 4], [5, 10], Strand.PLUS),
                3,
                3,
                0,
                [
                    SingleInterval(0, 3, Strand.PLUS),
                    CompoundInterval([3, 4], [5, 5], Strand.PLUS),
                    SingleInterval(5, 8, Strand.PLUS),
                ],
            ),
        ],
    )
    def test_scan_windows(self, location, window_size, step_size, start_pos, expected):
        assert list(location.scan_windows(window_size, step_size, start_pos)) == expected

    @pytest.mark.parametrize(
        "location,window_size,step_size,start_pos,expected_exception",
        [
            # Window size is greater than location length
            (CompoundInterval([0, 10], [5, 15], Strand.PLUS), 11, 1, 0, ValueError),
            # Location is unstranded
            (
                CompoundInterval([0, 10], [5, 15], Strand.UNSTRANDED),
                1,
                1,
                0,
                InvalidStrandException,
            ),
            # Window size is 0
            (CompoundInterval([0, 10], [5, 15], Strand.PLUS), 0, 1, 0, ValueError),
            # Step size is 0
            (CompoundInterval([0, 10], [5, 15], Strand.PLUS), 1, 0, 0, ValueError),
            # Start pos is negative
            (CompoundInterval([0, 10], [5, 15], Strand.PLUS), 1, 1, -1, ValueError),
            # Start pos is greater than location length
            (CompoundInterval([0, 10], [5, 15], Strand.PLUS), 1, 1, 20, ValueError),
            # Start pos + window size is greater than location length
            (CompoundInterval([0, 10], [5, 15], Strand.PLUS), 5, 1, 8, ValueError),
        ],
    )
    def test_scan_windows_error(self, location, window_size, step_size, start_pos, expected_exception):
        with pytest.raises(expected_exception):
            list(location.scan_windows(window_size, step_size, start_pos))

    @pytest.mark.parametrize(
        "location,other,match_strand,expected",
        [
            # Blocks interleaved, no overlap
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                CompoundInterval([5, 15], [10, 20], Strand.PLUS),
                True,
                False,
            ),
            # Different strands, do match strand
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                CompoundInterval([0, 10], [5, 15], Strand.UNSTRANDED),
                True,
                False,
            ),
            # Different strands, do not match strand
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                SingleInterval(0, 20, Strand.MINUS),
                False,
                True,
            ),
            # Middle blocks overlap, SingleInterval
            (
                CompoundInterval([0, 20, 40], [10, 30, 50], Strand.MINUS),
                SingleInterval(29, 32, Strand.MINUS),
                True,
                True,
            ),
            # Middle blocks overlap, CompoundInterval
            (
                CompoundInterval([0, 20, 40], [10, 30, 50], Strand.MINUS),
                CompoundInterval([11, 18, 70], [15, 22, 80], Strand.MINUS),
                True,
                True,
            ),
            # Ends overlap, SingleInterval
            (
                CompoundInterval([0, 20, 40], [10, 30, 50], Strand.MINUS),
                SingleInterval(49, 60, Strand.MINUS),
                True,
                True,
            ),
            # Ends overlap, CompoundInterval
            (
                CompoundInterval([0, 20, 40], [10, 30, 50], Strand.MINUS),
                CompoundInterval([49, 70], [55, 75], Strand.MINUS),
                True,
                True,
            ),
            # No overlap
            (
                CompoundInterval([0, 20, 40], [10, 30, 50], Strand.MINUS),
                CompoundInterval([50, 70], [55, 75], Strand.MINUS),
                True,
                False,
            ),
            # Same coordinates and strand, different parents
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS, parent="parent1"),
                CompoundInterval([0, 10], [5, 15], Strand.PLUS, parent="parent2"),
                False,
                False,
            ),
        ],
    )
    def test_has_overlap(self, location, other, match_strand, expected):
        assert location.has_overlap(other, match_strand) is expected
        assert other.has_overlap(location, match_strand) is expected

    @pytest.mark.parametrize(
        "location,expected",
        [
            # One block
            (
                CompoundInterval([5], [8], Strand.PLUS),
                SingleInterval(5, 8, Strand.PLUS),
            ),
            # Two adjacent blocks
            (
                CompoundInterval([5, 8], [8, 10], Strand.MINUS),
                SingleInterval(5, 10, Strand.MINUS),
            ),
            # Three adjacent blocks
            (
                CompoundInterval([5, 8, 10], [8, 10, 12], Strand.PLUS),
                SingleInterval(5, 12, Strand.PLUS),
            ),
            # Two sets of two adjacent blocks, one non-adjacent
            (
                CompoundInterval([5, 8, 15, 25, 30], [8, 10, 20, 30, 32], Strand.PLUS),
                CompoundInterval([5, 15, 25], [10, 20, 32], Strand.PLUS),
            ),
            # Has an empty block
            (
                CompoundInterval([0, 10, 20], [5, 10, 25], Strand.PLUS),
                CompoundInterval([0, 20], [5, 25], Strand.PLUS),
            ),
            # No adjacent blocks
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
            ),
            # Overlapping and adjacent blocks
            (
                CompoundInterval([5, 8, 10, 20], [9, 10, 12, 30], Strand.PLUS),
                CompoundInterval([5, 8, 20], [9, 12, 30], Strand.PLUS),
            ),
        ],
    )
    def test_optimize_blocks(self, location, expected):
        assert location.optimize_blocks() == expected

    @pytest.mark.parametrize(
        "location,expected",
        [
            # One block
            (
                CompoundInterval([5], [8], Strand.PLUS),
                [],
            ),
            # Two blocks
            (
                CompoundInterval([10, 20], [15, 25], Strand.PLUS),
                [SingleInterval(15, 20, Strand.PLUS)],
            ),
            # Three blocks, plus strand
            (
                CompoundInterval([10, 20, 30], [15, 25, 35], Strand.PLUS),
                [
                    SingleInterval(15, 20, Strand.PLUS),
                    SingleInterval(25, 30, Strand.PLUS),
                ],
            ),
            # Four blocks with parent, minus strand
            (
                CompoundInterval([10, 20, 30, 40], [15, 25, 35, 45], Strand.MINUS, parent="parent"),
                [
                    SingleInterval(35, 40, Strand.MINUS, parent="parent"),
                    SingleInterval(25, 30, Strand.MINUS, parent="parent"),
                    SingleInterval(15, 20, Strand.MINUS, parent="parent"),
                ],
            ),
            # Three adjacent blocks
            (
                CompoundInterval([5, 8, 10], [8, 10, 12], Strand.PLUS),
                [],
            ),
            # Two sets of two adjacent blocks, one non-adjacent
            (
                CompoundInterval([5, 8, 15, 25, 30], [8, 10, 20, 30, 32], Strand.MINUS),
                [
                    SingleInterval(20, 25, Strand.MINUS),
                    SingleInterval(10, 15, Strand.MINUS),
                ],
            ),
            # Has an empty block between two nonempty blocks
            (
                CompoundInterval([0, 10, 20], [5, 10, 25], Strand.PLUS),
                [SingleInterval(5, 20, Strand.PLUS)],
            ),
        ],
    )
    def test_gap_list(self, location, expected):
        assert location.gap_list() == expected

    @pytest.mark.parametrize(
        "location,expected",
        [
            # One block
            (
                CompoundInterval([5], [8], Strand.PLUS),
                EmptyLocation(),
            ),
            # Two blocks
            (
                CompoundInterval([10, 20], [15, 25], Strand.PLUS),
                CompoundInterval([15], [20], Strand.PLUS),
            ),
            # Four blocks with parent, minus strand
            (
                CompoundInterval([10, 20, 30, 40], [15, 25, 35, 45], Strand.MINUS, parent="parent"),
                CompoundInterval([15, 25, 35], [20, 30, 40], Strand.MINUS, parent="parent"),
            ),
            # Three adjacent blocks
            (
                CompoundInterval([5, 8, 10], [8, 10, 12], Strand.PLUS),
                EmptyLocation(),
            ),
            # Two sets of two adjacent blocks, one non-adjacent
            (
                CompoundInterval([5, 8, 15, 25, 30], [8, 10, 20, 30, 32], Strand.MINUS),
                CompoundInterval([10, 20], [15, 25], Strand.MINUS),
            ),
            # Has an empty block between two nonempty blocks
            (
                CompoundInterval([0, 10, 20], [5, 10, 25], Strand.PLUS),
                CompoundInterval([5], [20], Strand.PLUS),
            ),
        ],
    )
    def test_gaps_location(self, location, expected):
        assert location.gaps_location() == expected

    @pytest.mark.parametrize(
        "location,expected",
        [
            # Has parent with attributes; minus strand
            (
                CompoundInterval(
                    [10, 20],
                    [17, 22],
                    Strand.MINUS,
                    parent=Parent(
                        id="parent",
                        location=CompoundInterval([10, 20], [17, 22], Strand.MINUS),
                    ),
                ),
                CompoundInterval([10, 15], [12, 22], Strand.PLUS, parent="parent"),
            ),
            # Plus strand
            (
                CompoundInterval([0, 20, 40], [10, 22, 45], Strand.PLUS),
                CompoundInterval([0, 23, 35], [5, 25, 45], Strand.MINUS),
            ),
            # One block, plus strand
            (
                CompoundInterval([5], [10], Strand.PLUS),
                CompoundInterval([5], [10], Strand.MINUS),
            ),
            # Same blocks both orientations
            (
                CompoundInterval([0, 10], [5, 15], Strand.MINUS),
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
            ),
            # Unstranded, multiple blocks
            (
                CompoundInterval([0, 20, 40], [10, 22, 45], Strand.UNSTRANDED),
                CompoundInterval([0, 23, 35], [5, 25, 45], Strand.UNSTRANDED),
            ),
        ],
    )
    def test_reverse(self, location, expected):
        assert location.reverse() == expected

    @pytest.mark.parametrize(
        "location,expected",
        [
            # Has parent with attributes; minus strand
            (
                CompoundInterval(
                    [10, 20],
                    [17, 22],
                    Strand.MINUS,
                    parent=Parent(
                        id="parent",
                        location=CompoundInterval([10, 20], [17, 22], Strand.MINUS),
                    ),
                ),
                CompoundInterval([10, 20], [17, 22], Strand.PLUS, parent="parent"),
            ),
            # Plus strand
            (
                CompoundInterval([0, 20, 40], [10, 22, 45], Strand.PLUS),
                CompoundInterval([0, 20, 40], [10, 22, 45], Strand.MINUS),
            ),
            # Unstranded
            (
                CompoundInterval([0, 20, 40], [10, 22, 45], Strand.UNSTRANDED),
                CompoundInterval([0, 20, 40], [10, 22, 45], Strand.UNSTRANDED),
            ),
        ],
    )
    def test_reverse_strand(self, location, expected):
        assert location.reverse_strand() == expected

    def test_reset_strand(self):
        # Has parent with attributes
        location = CompoundInterval(
            [10, 20],
            [17, 22],
            Strand.MINUS,
            parent=Parent(
                id="parent",
                location=CompoundInterval([10, 20], [17, 22], Strand.MINUS),
                sequence_type="chr",
            ),
        )
        assert location.reset_strand(Strand.PLUS) == CompoundInterval(
            [10, 20],
            [17, 22],
            Strand.PLUS,
            parent=Parent(id="parent", sequence_type="chr"),
        )

    @pytest.mark.parametrize(
        "location,new_parent,expected",
        [
            # Old and new parents have all attributes
            (
                CompoundInterval(
                    [0, 3],
                    [1, 4],
                    Strand.PLUS,
                    parent=Parent(id="old_parent", sequence=Sequence("AAAA", Alphabet.NT_STRICT)),
                ),
                Parent(id="new_parent", sequence=Sequence("TTTTT", Alphabet.NT_STRICT)),
                CompoundInterval(
                    [0, 3],
                    [1, 4],
                    Strand.PLUS,
                    parent=Parent(id="new_parent", sequence=Sequence("TTTTT", Alphabet.NT_STRICT)),
                ),
            ),
            # No original parent
            (
                CompoundInterval([0, 3], [1, 4], Strand.PLUS),
                Parent(id="new_parent", sequence=Sequence("TTTTT", Alphabet.NT_STRICT)),
                CompoundInterval(
                    [0, 3],
                    [1, 4],
                    Strand.PLUS,
                    parent=Parent(id="new_parent", sequence=Sequence("TTTTT", Alphabet.NT_STRICT)),
                ),
            ),
            # No new parent
            (
                CompoundInterval(
                    [0, 3],
                    [1, 4],
                    Strand.PLUS,
                    parent=Parent(id="old_parent", sequence=Sequence("AAAA", Alphabet.NT_STRICT)),
                ),
                None,
                CompoundInterval([0, 3], [1, 4], Strand.PLUS),
            ),
        ],
    )
    def test_reset_parent(self, location, new_parent, expected):
        assert location.reset_parent(new_parent) == expected

    def test_reset_parent_error(self):
        # New parent is shorter than location
        location = CompoundInterval(
            [0, 3],
            [1, 4],
            Strand.PLUS,
            parent=Parent(id="old_parent", sequence=Sequence("AAAA", Alphabet.NT_STRICT)),
        )
        with pytest.raises(InvalidPositionException):
            location.reset_parent(Parent(id="new_parent", sequence=Sequence("T", Alphabet.NT_STRICT)))

    @pytest.mark.parametrize(
        "location,expected",
        [
            # Empty block
            (CompoundInterval([6], [6], Strand.PLUS), 1),
            # One block
            (CompoundInterval([0], [3], Strand.PLUS), 1),
            # Multiple blocks
            (CompoundInterval([0, 5], [3, 7], Strand.PLUS), 2),
            # Multiple blocks, some adjacent to each other
            (CompoundInterval([0, 5, 10], [5, 8, 12], Strand.PLUS), 3),
        ],
    )
    def test_num_blocks(self, location, expected):
        assert location.num_blocks == expected

    @pytest.mark.parametrize(
        "location,other,distance_type,expected",
        [
            # Both have parents, SingleInterval
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS, parent="parent"),
                SingleInterval(20, 30, Strand.MINUS, parent="parent"),
                DistanceType.ENDS,
                15,
            ),
            # Overlapping, different strands, inner distance
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                CompoundInterval([13, 20], [16, 25], Strand.UNSTRANDED),
                DistanceType.INNER,
                0,
            ),
            # Overlapping, outer distance
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                CompoundInterval([13, 20], [16, 25], Strand.UNSTRANDED),
                DistanceType.OUTER,
                25,
            ),
            # Disjoint, inner distance, SingleInterval
            (
                CompoundInterval([5, 10], [7, 15], Strand.PLUS),
                SingleInterval(0, 3, Strand.PLUS),
                DistanceType.INNER,
                2,
            ),
            # Disjoint, adjacent, inner distance
            (
                CompoundInterval([5, 10], [7, 15], Strand.PLUS),
                CompoundInterval([15, 20], [17, 25], Strand.PLUS),
                DistanceType.INNER,
                0,
            ),
            # Disjoint, adjacent, outer distance
            (
                CompoundInterval([5, 10], [7, 15], Strand.PLUS),
                CompoundInterval([15, 20], [17, 25], Strand.PLUS),
                DistanceType.OUTER,
                20,
            ),
            # Disjoint, starts distance
            (
                CompoundInterval([5, 10], [7, 15], Strand.PLUS),
                CompoundInterval([15, 20], [17, 25], Strand.PLUS),
                DistanceType.STARTS,
                10,
            ),
            # Overlapping, ends distance
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                CompoundInterval([13, 20], [16, 25], Strand.UNSTRANDED),
                DistanceType.ENDS,
                10,
            ),
            # Interleaved blocks no overlap, outer distance
            (
                CompoundInterval([0, 20, 40], [10, 30, 50], Strand.PLUS),
                CompoundInterval([15, 35], [17, 39], Strand.PLUS),
                DistanceType.OUTER,
                39,
            ),
            # Interleaved blocks no overlap, inner distance
            (
                CompoundInterval([0, 20, 40], [10, 30, 50], Strand.PLUS),
                CompoundInterval([15, 35], [17, 39], Strand.PLUS),
                DistanceType.INNER,
                1,
            ),
            # Same coordinates, outer distance
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                DistanceType.OUTER,
                15,
            ),
        ],
    )
    def test_distance_to(self, location, other, distance_type, expected):
        assert location.distance_to(other, distance_type) == expected

    # TODO: RE-enable this test if we decide to throw a ValueError for this kind of comparison
    #    for compound intervals (and presumably for simple intervals as well!)
    # def test_distance_to_error(self):
    #     # Different parents
    #     with pytest.raises(ValueError):
    #         CompoundInterval([0, 10], [5, 15], Strand.PLUS, parent="parent1").distance_to(
    #             CompoundInterval([0, 10], [5, 15], Strand.PLUS, parent="parent2")
    #         )

    @pytest.mark.parametrize(
        "location,sequence_type,expected",
        [
            # Parent
            (
                CompoundInterval(
                    [0],
                    [10],
                    Strand.PLUS,
                    parent=Parent(id="parent", sequence_type="chr"),
                ),
                "chr",
                Parent(
                    id="parent",
                    sequence_type="chr",
                    location=CompoundInterval([0], [10], Strand.PLUS),
                ),
            ),
            # Grandparent
            (
                CompoundInterval(
                    [0],
                    [10],
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
            # Grandparent; parent has no type
            (
                CompoundInterval(
                    [0],
                    [10],
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
    def test_first_ancestor_of_type(self, location, sequence_type, expected):
        assert location.first_ancestor_of_type(sequence_type) == expected

    @pytest.mark.parametrize(
        "location,sequence_type",
        [
            # No parent
            (CompoundInterval([0], [10], Strand.PLUS), "seq_chunk"),
            # Parent has no type
            (
                CompoundInterval([0], [10], Strand.PLUS, parent="parent"),
                "seq_chunk",
            ),
            # Parent has wrong type, no additional ancestors
            (
                CompoundInterval(
                    [0],
                    [10],
                    Strand.PLUS,
                    parent=Parent(id="parent", sequence_type="chr"),
                ),
                "seq_chunk",
            ),
            # Two ancestors with types, no matching type
            (
                CompoundInterval(
                    [0],
                    [10],
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="chr",
                        parent=Parent(
                            id="grandparent",
                            sequence_type="seqtype",
                        ),
                    ),
                ),
                "seq_chunk",
            ),
        ],
    )
    def test_first_ancestor_of_type_error(self, location, sequence_type):
        with pytest.raises(NoSuchAncestorException):
            location.first_ancestor_of_type(sequence_type)

    @pytest.mark.parametrize(
        "location,sequence_type,expected",
        [
            # Parent
            (
                CompoundInterval(
                    [0],
                    [10],
                    Strand.PLUS,
                    parent=Parent(id="parent", sequence_type="chr"),
                ),
                "chr",
                True,
            ),
            # Grandparent
            (
                CompoundInterval(
                    [0],
                    [10],
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
            # Grandparent; parent has no type
            (
                CompoundInterval(
                    [0],
                    [10],
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
            (CompoundInterval([0], [10], Strand.PLUS), "seq_chunk", False),
            # Parent has no type
            (
                CompoundInterval([0], [10], Strand.PLUS, parent="parent"),
                "seq_chunk",
                False,
            ),
            # Parent has wrong type, no additional ancestors
            (
                CompoundInterval(
                    [0],
                    [10],
                    Strand.PLUS,
                    parent=Parent(id="parent", sequence_type="chr"),
                ),
                "seq_chunk",
                False,
            ),
            # Two ancestors with types, no matching type
            (
                CompoundInterval(
                    [0],
                    [10],
                    Strand.PLUS,
                    parent=Parent(
                        id="parent",
                        sequence_type="chr",
                        parent=Parent(
                            id="grandparent",
                            sequence_type="seqtype",
                        ),
                    ),
                ),
                "seq_chunk",
                False,
            ),
        ],
    )
    def test_has_ancestor_of_type(self, location, sequence_type, expected):
        assert location.has_ancestor_of_type(sequence_type) is expected

    @pytest.mark.parametrize(
        "location,sequence_type,expected",
        [
            # Parent
            (
                CompoundInterval(
                    [5],
                    [10],
                    Strand.MINUS,
                    parent=Parent(sequence_type="chr"),
                ),
                "chr",
                CompoundInterval(
                    [5],
                    [10],
                    Strand.MINUS,
                    parent=Parent(sequence_type="chr"),
                ),
            ),
            # Grandparent, parent has no type
            (
                CompoundInterval(
                    [5],
                    [10],
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
                CompoundInterval(
                    [2],
                    [4],
                    Strand.PLUS,
                    Parent(
                        id="parent",
                        sequence_type="chr",
                        strand=Strand.PLUS,
                        location=CompoundInterval([2], [4], Strand.PLUS),
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
                CompoundInterval(
                    [2],
                    [4],
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
                CompoundInterval(
                    [2],
                    [4],
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
                CompoundInterval(
                    [2],
                    [4],
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
                CompoundInterval(
                    [2],
                    [4],
                    Strand.MINUS,
                    Parent(
                        id="parent",
                        sequence_type="seq_type",
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
        ],
    )
    def test_lift_over_to_first_ancestor_of_type(self, location, sequence_type, expected):
        assert location.lift_over_to_first_ancestor_of_type(sequence_type) == expected

    @pytest.mark.parametrize(
        "location,sequence_type,expected_exception",
        [
            # No parent
            (
                CompoundInterval([5], [10], Strand.PLUS),
                "seq_chunk",
                NoSuchAncestorException,
            ),
            # Parent has no type
            (
                CompoundInterval(
                    [5],
                    [10],
                    Strand.PLUS,
                    parent=Parent(sequence_type="seq_chunk"),
                ),
                "unknown",
                NoSuchAncestorException,
            ),
            # Parent has wrong type, no additional ancestors
            (
                CompoundInterval(
                    [5],
                    [10],
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
                CompoundInterval(
                    [2],
                    [4],
                    Strand.MINUS,
                    Parent(
                        id="parent",
                        sequence_type="seq_chunk",
                        parent=Parent(id="grandparent", sequence_type="chr"),
                    ),
                ),
                "chr",
                ValueError,
            ),
            # Parent has no location on grandparent
            (
                CompoundInterval(
                    [2],
                    [4],
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
                ValueError,
            ),
        ],
    )
    def test_lift_over_to_first_ancestor_of_type_error(self, location, sequence_type, expected_exception):
        with pytest.raises(expected_exception):
            location.lift_over_to_first_ancestor_of_type(sequence_type)

    @pytest.mark.parametrize(
        "location,shift,expected",
        [
            # Positive shift
            (
                CompoundInterval([10, 30], [11, 35], Strand.MINUS),
                5,
                CompoundInterval([15, 35], [16, 40], Strand.MINUS),
            ),
            # Negative shift
            (
                CompoundInterval([10, 30], [11, 35], Strand.MINUS),
                -5,
                CompoundInterval([5, 25], [6, 30], Strand.MINUS),
            ),
            # Has parent
            (
                CompoundInterval(
                    [10, 30],
                    [11, 35],
                    Strand.MINUS,
                    parent=Parent(
                        id="parent",
                        location=CompoundInterval([10, 30], [11, 35], Strand.MINUS),
                    ),
                ),
                2,
                CompoundInterval([12, 32], [13, 37], Strand.MINUS, parent="parent"),
            ),
        ],
    )
    def test_shift_position(self, location, shift, expected):
        assert location.shift_position(shift) == expected

    @pytest.mark.parametrize(
        "location,shift,expected_exception",
        [
            # Off beginning of sequence
            (
                CompoundInterval([10, 30], [11, 35], Strand.MINUS),
                -20,
                InvalidPositionException,
            ),
            # Off end of sequence
            (
                CompoundInterval(
                    [0, 3],
                    [2, 5],
                    Strand.MINUS,
                    parent=Sequence("AAAAAA", Alphabet.NT_STRICT),
                ),
                10,
                InvalidPositionException,
            ),
        ],
    )
    def test_shift_position_error(self, location, shift, expected_exception):
        with pytest.raises(expected_exception):
            location.shift_position(shift)

    @pytest.mark.parametrize(
        "location,sequence,expected",
        [
            # No parent
            (
                CompoundInterval([0], [1], Strand.PLUS),
                Sequence("AA", Alphabet.NT_STRICT),
                False,
            ),
            # Has parent with given sequence
            (
                CompoundInterval(
                    [0],
                    [1],
                    Strand.PLUS,
                    parent=Sequence("AA", Alphabet.NT_STRICT),
                ),
                Sequence("AA", Alphabet.NT_STRICT),
                True,
            ),
            # Has parent with sequence with one different attribute
            (
                CompoundInterval(
                    [0],
                    [1],
                    Strand.PLUS,
                    parent=Sequence("AA", Alphabet.NT_STRICT),
                ),
                Sequence("AA", Alphabet.NT_STRICT, id="id"),
                False,
            ),
            # Grandparent
            (
                CompoundInterval(
                    [0],
                    [1],
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
                CompoundInterval(
                    [0],
                    [1],
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
    def test_has_ancestor_sequence(self, location, sequence, expected):
        assert location.has_ancestor_sequence(sequence) is expected

    @pytest.mark.parametrize(
        "location,sequence,expected",
        [
            # Parent
            (
                CompoundInterval(
                    [1],
                    [2],
                    Strand.PLUS,
                    parent=Sequence("AAA", Alphabet.NT_STRICT),
                ),
                Sequence("AAA", Alphabet.NT_STRICT),
                CompoundInterval(
                    [1],
                    [2],
                    Strand.PLUS,
                    parent=Sequence("AAA", Alphabet.NT_STRICT),
                ),
            ),
            # Grandparent, plus/minus
            (
                CompoundInterval(
                    [1],
                    [2],
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
            # Grandparent, parent has no sequence
            (
                CompoundInterval(
                    [1],
                    [2],
                    Strand.PLUS,
                    parent=Parent(
                        parent=Parent(
                            location=SingleInterval(2, 5, Strand.MINUS),
                            sequence=Sequence("TTTTT", Alphabet.NT_STRICT),
                        )
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
    def test_lift_over_to_sequence(self, location, sequence, expected):
        assert location.lift_over_to_sequence(sequence) == expected

    @pytest.mark.parametrize(
        "location,sequence,expected_exception",
        [
            # Location not contiguous
            (
                CompoundInterval(
                    [1, 5],
                    [2, 7],
                    Strand.PLUS,
                    parent=Sequence("AAAAAAAA", Alphabet.NT_STRICT),
                ),
                Sequence("AAA", Alphabet.NT_STRICT),
                ValueError,
            ),
            # No parent
            (
                CompoundInterval([0], [1], Strand.PLUS),
                Sequence("AA", Alphabet.NT_STRICT),
                NoSuchAncestorException,
            ),
            # Has parent with sequence with different id
            (
                CompoundInterval(
                    [0],
                    [1],
                    Strand.PLUS,
                    parent=Sequence("AA", Alphabet.NT_STRICT),
                ),
                Sequence("AA", Alphabet.NT_STRICT, id="id"),
                NoSuchAncestorException,
            ),
            # Parent doesn't match, parent has empty parent
            (
                CompoundInterval(
                    [0],
                    [1],
                    Strand.PLUS,
                    parent=Parent(sequence=Sequence("AA", Alphabet.NT_STRICT), parent=Parent()),
                ),
                Sequence("AA", Alphabet.NT_STRICT, id="id"),
                NoSuchAncestorException,
            ),
        ],
    )
    def test_lift_over_to_sequence_error(self, location, sequence, expected_exception):
        with pytest.raises(expected_exception):
            location.lift_over_to_sequence(sequence)

    @pytest.mark.parametrize(
        "single_interval,compound_interval,expected",
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
            # Single interval overlaps 3 blocks
            (
                SingleInterval(15, 60, Strand.PLUS),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.PLUS),
                SingleInterval(10, 60, Strand.PLUS),
            ),
            # Single interval contains entire compound interval
            (
                SingleInterval(0, 80, Strand.UNSTRANDED),
                CompoundInterval([10, 30, 50], [20, 40, 60], Strand.UNSTRANDED),
                SingleInterval(0, 80, Strand.UNSTRANDED),
            ),
        ],
    )
    def test_union_single_interval(self, single_interval, compound_interval, expected):
        assert compound_interval.union(single_interval) == expected

    @pytest.mark.parametrize(
        "compound_interval,single_interval,expected_exception",
        [
            # Both empty
            (
                CompoundInterval([5], [5], Strand.PLUS),
                SingleInterval(10, 10, Strand.PLUS),
                ValueError,
            ),
            # Opposite strands
            (
                CompoundInterval([5], [10], Strand.PLUS),
                SingleInterval(10, 15, Strand.UNSTRANDED),
                ValueError,
            ),
            # Different parents
            (
                CompoundInterval([0], [5], Strand.PLUS, parent="parent1"),
                SingleInterval(3, 10, Strand.PLUS, parent="parent2"),
                ValueError,
            ),
            # One has no parent
            (
                CompoundInterval([0], [5], Strand.PLUS, parent="parent1"),
                SingleInterval(3, 10, Strand.PLUS),
                ValueError,
            ),
        ],
    )
    def test_union_single_interval_error(self, compound_interval, single_interval, expected_exception):
        with pytest.raises(expected_exception):
            compound_interval.union(single_interval)

    @pytest.mark.parametrize(
        "location1,location2,expected",
        [
            # Location 1 empty
            (
                CompoundInterval([5], [5], Strand.PLUS),
                CompoundInterval([3], [4], Strand.PLUS),
                SingleInterval(3, 4, Strand.PLUS),
            ),
            # Location 2 empty
            (
                CompoundInterval([10, 20], [15, 30], Strand.MINUS),
                CompoundInterval([7], [7], Strand.MINUS),
                CompoundInterval([10, 20], [15, 30], Strand.MINUS),
            ),
            # Disjoint, not touching, both have parents
            (
                CompoundInterval([10, 20], [15, 30], Strand.MINUS, parent="parent"),
                CompoundInterval([0], [5], Strand.MINUS, parent="parent"),
                CompoundInterval([0, 10, 20], [5, 15, 30], Strand.MINUS, parent="parent"),
            ),
            # Disjoint, adjacent
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                CompoundInterval([15, 25], [18, 30], Strand.PLUS),
                CompoundInterval([0, 10, 25], [5, 18, 30], Strand.PLUS),
            ),
            # Blocks perfectly interleaved
            (
                CompoundInterval([0, 10, 20], [5, 15, 25], Strand.UNSTRANDED),
                CompoundInterval([5, 15], [10, 20], Strand.UNSTRANDED),
                SingleInterval(0, 25, Strand.UNSTRANDED),
            ),
            # Blocks interleaved, not touching
            (
                CompoundInterval([0, 10, 20], [5, 15, 25], Strand.UNSTRANDED),
                CompoundInterval([6, 16], [9, 19], Strand.UNSTRANDED),
                CompoundInterval([0, 6, 10, 16, 20], [5, 9, 15, 19, 25], Strand.UNSTRANDED),
            ),
            # Several blocks overlap on different sides
            (
                CompoundInterval([0, 10, 20], [5, 15, 25], Strand.PLUS),
                CompoundInterval([0, 3, 18], [1, 7, 26], Strand.PLUS),
                CompoundInterval([0, 10, 18], [7, 15, 26], Strand.PLUS),
            ),
            # Location 1 contained in a block of location 2
            (
                CompoundInterval([5, 10], [8, 13], Strand.MINUS),
                CompoundInterval([3, 20], [15, 25], Strand.MINUS),
                CompoundInterval([3, 20], [15, 25], Strand.MINUS),
            ),
            # Location 2 contained in an intron of location 1
            (
                CompoundInterval([0, 10, 30], [5, 15, 35], Strand.PLUS),
                CompoundInterval([20, 25], [22, 30], Strand.PLUS),
                CompoundInterval([0, 10, 20, 25], [5, 15, 22, 35], Strand.PLUS),
            ),
            # identical
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
            ),
            # enclosed
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                CompoundInterval([0], [20], Strand.PLUS),
                SingleInterval(0, 20, Strand.PLUS),
            ),
            # overlapping
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                CompoundInterval([3, 10], [8, 20], Strand.PLUS),
                CompoundInterval([0, 10], [8, 20], Strand.PLUS),
            ),
            # disjoint
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                CompoundInterval([30], [40], Strand.PLUS),
                CompoundInterval([0, 10, 30], [5, 20, 40], Strand.PLUS),
            ),
        ],
    )
    def test_union_compound_interval(self, location1, location2, expected):
        assert location1.union(location2) == expected

    @pytest.mark.parametrize(
        "interval,expected",
        [  # this is a real example of a -1 frameshift gene in E. coli
            (
                CompoundInterval(
                    [4507196, 4509008, 4509478, 4509795],
                    [4507451, 4509479, 4509793, 4509800],
                    Strand.PLUS,
                ),
                CompoundInterval(
                    [4507196, 4509008, 4509795],
                    [4507451, 4509793, 4509800],
                    Strand.PLUS,
                ),
            ),
            # no overlapping blocks comes back unchanged
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
            ),
        ],
    )
    def test_merge_overlapping_interval(self, interval, expected):
        assert interval.merge_overlapping() == expected

    @pytest.mark.parametrize(
        "location1,location2,expected_exception",
        [
            # Opposite strands
            (
                CompoundInterval([5], [6], Strand.MINUS),
                CompoundInterval([3], [5], Strand.PLUS),
                ValueError,
            ),
            # Different parents
            (
                CompoundInterval([0], [5], Strand.PLUS, parent="parent1"),
                CompoundInterval([3], [7], Strand.PLUS, parent="parent2"),
                ValueError,
            ),
            # One has no parent
            (
                CompoundInterval([0], [5], Strand.PLUS, parent="parent1"),
                CompoundInterval([3], [7], Strand.PLUS),
                ValueError,
            ),
        ],
    )
    def test_union_compound_interval_error(self, location1, location2, expected_exception):
        with pytest.raises(expected_exception):
            location1.union(location2)

    @pytest.mark.parametrize(
        "location1,location2,match_strand,expected",
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
                CompoundInterval([5], [10], Strand.PLUS),
                SingleInterval(5, 10, Strand.MINUS),
                False,
                SingleInterval(5, 10, Strand.PLUS),
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
    def test_intersection_single_interval(self, location1, location2, match_strand, expected):
        assert location1.intersection(location2, match_strand) == expected

    @pytest.mark.parametrize(
        "location1,location2,match_strand,expected",
        [
            # Different strands
            (
                CompoundInterval([5], [10], Strand.PLUS),
                CompoundInterval([5], [10], Strand.MINUS),
                True,
                EmptyLocation(),
            ),
            # Different strands, match_strand=False
            (
                CompoundInterval([5], [10], Strand.PLUS),
                CompoundInterval([5], [10], Strand.MINUS),
                False,
                SingleInterval(5, 10, Strand.PLUS),
            ),
            # No overlap
            (
                CompoundInterval([5], [10], Strand.PLUS),
                CompoundInterval([10], [20], Strand.PLUS),
                True,
                EmptyLocation(),
            ),
            # First contains second
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                CompoundInterval([3, 12], [4, 16], Strand.PLUS),
                True,
                CompoundInterval([3, 12], [4, 16], Strand.PLUS),
            ),
            # Second contains first
            (
                CompoundInterval([3, 12], [4, 16], Strand.PLUS),
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                True,
                CompoundInterval([3, 12], [4, 16], Strand.PLUS),
            ),
            # Multiple blocks overlap
            (
                CompoundInterval([0, 10, 20], [5, 15, 25], Strand.MINUS),
                CompoundInterval([3, 17], [12, 30], Strand.MINUS),
                True,
                CompoundInterval([3, 10, 20], [5, 12, 25], Strand.MINUS),
            ),
            # Both have parents and locations
            (
                CompoundInterval(
                    [5],
                    [10],
                    Strand.PLUS,
                    parent=Parent(id="parent", location=CompoundInterval([5], [10], Strand.PLUS)),
                ),
                CompoundInterval(
                    [5],
                    [8],
                    Strand.PLUS,
                    parent=Parent(id="parent", location=CompoundInterval([5], [8], Strand.PLUS)),
                ),
                True,
                SingleInterval(5, 8, Strand.PLUS, parent="parent"),
            ),
            # Different parents
            (
                CompoundInterval([5], [8], Strand.PLUS, parent="parent1"),
                CompoundInterval([5], [8], Strand.PLUS, parent="parent2"),
                True,
                EmptyLocation(),
            ),
        ],
    )
    def test_intersection_compound_interval(self, location1, location2, match_strand, expected):
        assert location1.intersection(location2, match_strand) == expected

    @pytest.mark.parametrize(
        "location1,location2,match_strand,expected",
        [
            # No overlap
            (
                CompoundInterval([5], [10], Strand.PLUS),
                SingleInterval(10, 20, Strand.PLUS),
                True,
                SingleInterval(5, 10, Strand.PLUS),
            ),
            # Different strands
            (
                CompoundInterval([5], [10], Strand.PLUS),
                SingleInterval(5, 10, Strand.MINUS),
                True,
                SingleInterval(5, 10, Strand.PLUS),
            ),
            # Different strands, match_strand=False
            (
                CompoundInterval([5], [10], Strand.PLUS),
                SingleInterval(5, 8, Strand.MINUS),
                False,
                SingleInterval(8, 10, Strand.PLUS),
            ),
            # Different parents
            (
                CompoundInterval([5], [10], Strand.PLUS, parent="parent1"),
                SingleInterval(5, 10, Strand.PLUS, parent="parent2"),
                True,
                SingleInterval(5, 10, Strand.PLUS, parent="parent1"),
            ),
            # Single interval contains compound interval
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                SingleInterval(0, 20, Strand.PLUS),
                True,
                EmptyLocation(),
            ),
            # Compound interval contains single interval
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                SingleInterval(3, 4, Strand.PLUS),
                True,
                CompoundInterval([0, 4, 10], [3, 5, 15], Strand.PLUS),
            ),
            # Multiple blocks overlap
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                SingleInterval(3, 12, Strand.PLUS),
                True,
                CompoundInterval([0, 12], [3, 15], Strand.PLUS),
            ),
            # Single interval contains one block of compound interval
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                SingleInterval(3, 20, Strand.PLUS),
                True,
                SingleInterval(0, 3, Strand.PLUS),
            ),
        ],
    )
    def test_minus_single_interval(self, location1, location2, match_strand, expected):
        assert location1.minus(location2, match_strand) == expected

    @pytest.mark.parametrize(
        "location1,location2,match_strand,expected",
        [
            # No overlap
            (
                CompoundInterval([0], [5], Strand.PLUS),
                CompoundInterval([5], [10], Strand.PLUS),
                True,
                SingleInterval(0, 5, Strand.PLUS),
            ),
            # Different strands
            (
                CompoundInterval([0], [5], Strand.PLUS),
                CompoundInterval([0], [5], Strand.MINUS),
                True,
                SingleInterval(0, 5, Strand.PLUS),
            ),
            # Different strands, match_strand=False
            (
                CompoundInterval([0], [5], Strand.PLUS),
                CompoundInterval([0], [3], Strand.MINUS),
                False,
                SingleInterval(3, 5, Strand.PLUS),
            ),
            # Different parents
            (
                CompoundInterval(
                    [0],
                    [5],
                    Strand.PLUS,
                    parent=Parent(id="parent1", location=CompoundInterval([0], [5], Strand.PLUS)),
                ),
                CompoundInterval([0], [5], Strand.PLUS, parent="parent2"),
                True,
                SingleInterval(0, 5, Strand.PLUS, parent="parent1"),
            ),
            # First contains second
            (
                CompoundInterval([0, 20], [10, 30], Strand.PLUS),
                CompoundInterval([5, 25], [8, 27], Strand.PLUS),
                True,
                CompoundInterval([0, 8, 20, 27], [5, 10, 25, 30], Strand.PLUS),
            ),
            (
                CompoundInterval([0, 20], [10, 30], Strand.PLUS),
                CompoundInterval([0, 25], [10, 27], Strand.PLUS),
                True,
                CompoundInterval([20, 27], [25, 30], Strand.PLUS),
            ),
            # Second contains first
            (
                CompoundInterval([5, 25], [8, 27], Strand.PLUS),
                CompoundInterval([0, 20], [10, 30], Strand.PLUS),
                True,
                EmptyLocation(),
            ),
            (
                CompoundInterval([5, 20], [8, 30], Strand.PLUS),
                CompoundInterval([0, 20], [10, 30], Strand.PLUS),
                True,
                EmptyLocation(),
            ),
            # Multiple blocks overlap
            (
                CompoundInterval([0, 20], [10, 30], Strand.PLUS),
                CompoundInterval([8, 19], [15, 22], Strand.PLUS),
                True,
                CompoundInterval([0, 22], [8, 30], Strand.PLUS),
            ),
            # Second contains one full block of first
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                CompoundInterval([10], [15], Strand.PLUS),
                True,
                SingleInterval(0, 5, Strand.PLUS),
            ),
        ],
    )
    def test_minus_compound_interval(self, location1, location2, match_strand, expected):
        assert location1.minus(location2, match_strand) == expected

    @pytest.mark.parametrize(
        "location,extend_left,extend_right,expected",
        [
            # Plus strand one block
            (CompoundInterval([5], [10], Strand.PLUS), 2, 3, SingleInterval(3, 13, Strand.PLUS)),
            # Minus strand two blocks
            (
                CompoundInterval([5, 15], [10, 20], Strand.MINUS),
                2,
                3,
                CompoundInterval([3, 15], [10, 23], Strand.MINUS),
            ),
            # Unstranded
            (
                CompoundInterval([5, 15], [10, 20], Strand.UNSTRANDED),
                2,
                3,
                CompoundInterval([3, 15], [10, 23], Strand.UNSTRANDED),
            ),
            # Has parent
            (
                CompoundInterval([5], [10], Strand.PLUS, parent="parent"),
                2,
                3,
                SingleInterval(
                    3, 13, Strand.PLUS, parent=Parent(id="parent", location=SingleInterval(3, 13, Strand.PLUS))
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
            (CompoundInterval([5], [10], Strand.PLUS), -1, 1, ValueError),
            # extend_right negative
            (CompoundInterval([5], [10], Strand.PLUS), 1, -1, ValueError),
            # Extends left beyond parent
            (CompoundInterval([5], [10], Strand.PLUS), 6, 1, InvalidPositionException),
            # Extends right beyond parent
            (
                CompoundInterval([0], [5], Strand.PLUS, parent=Sequence("AAAAAA", Alphabet.NT_STRICT)),
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
            # Plus strand one block
            (CompoundInterval([5], [10], Strand.PLUS), 2, 3, SingleInterval(3, 13, Strand.PLUS)),
            # Minus strand two blocks
            (
                CompoundInterval([5, 15], [10, 20], Strand.MINUS),
                2,
                3,
                CompoundInterval([2, 15], [10, 22], Strand.MINUS),
            ),
            # Has parent
            (
                CompoundInterval([5], [10], Strand.PLUS, parent="parent"),
                2,
                3,
                SingleInterval(
                    3, 13, Strand.PLUS, parent=Parent(id="parent", location=SingleInterval(3, 13, Strand.PLUS))
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
            (CompoundInterval([5], [10], Strand.UNSTRANDED), 1, 1, InvalidStrandException),
            # extend_left negative
            (CompoundInterval([5], [10], Strand.PLUS), -1, 1, ValueError),
            # extend_right negative
            (CompoundInterval([5], [10], Strand.PLUS), 1, -1, ValueError),
            # Extends left beyond parent
            (CompoundInterval([5], [10], Strand.MINUS), 0, 6, InvalidPositionException),
            (CompoundInterval([5], [10], Strand.PLUS), 6, 0, InvalidPositionException),
            # Extends right beyond parent
            (
                CompoundInterval([0], [5], Strand.PLUS, parent=Sequence("AAAAAA", Alphabet.NT_STRICT)),
                0,
                10,
                InvalidPositionException,
            ),
            (
                CompoundInterval([0], [5], Strand.MINUS, parent=Sequence("AAAAAA", Alphabet.NT_STRICT)),
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
        "location1,location2,match_strand,expected",
        [
            # Contains, different strand, match strand = true
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                SingleInterval(0, 3, Strand.MINUS),
                True,
                False,
            ),
            # Contains, different strand, match strand = false
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                SingleInterval(0, 3, Strand.MINUS),
                False,
                True,
            ),
            # Contains, same strand, match strand = true
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                SingleInterval(0, 3, Strand.PLUS),
                True,
                True,
            ),
            # Overlaps, doesn't contain
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                SingleInterval(0, 6, Strand.PLUS),
                False,
                False,
            ),
            # No overlap
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                SingleInterval(20, 25, Strand.PLUS),
                False,
                False,
            ),
            # Same coordinates same parent
            (
                CompoundInterval([0], [5], Strand.PLUS, parent="parent"),
                SingleInterval(0, 5, Strand.PLUS, parent="parent"),
                False,
                True,
            ),
            # Same coordinates different parent
            (
                CompoundInterval([0], [5], Strand.PLUS, parent="parent1"),
                SingleInterval(0, 5, Strand.PLUS, parent="parent2"),
                False,
                False,
            ),
        ],
    )
    def test_contains_single_interval(self, location1, location2, match_strand, expected):
        assert location1.contains(location2, match_strand) is expected

    @pytest.mark.parametrize(
        "location1,location2,match_strand,expected",
        [
            # Contains, different strand, match strand = true
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                CompoundInterval([0, 10], [5, 20], Strand.MINUS),
                True,
                False,
            ),
            # Contains, different strand, match strand = false
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                CompoundInterval([0, 10], [5, 20], Strand.MINUS),
                False,
                True,
            ),
            # Contains, same strand, match strand = true
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                True,
                True,
            ),
            # Overlaps, doesn't contain
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                CompoundInterval([0, 10], [5, 21], Strand.PLUS),
                False,
                False,
            ),
            # No overlap
            (
                CompoundInterval([0, 10], [5, 20], Strand.PLUS),
                CompoundInterval([20], [30], Strand.PLUS),
                False,
                False,
            ),
            # Same coordinates same parent
            (
                CompoundInterval([0], [5], Strand.PLUS, parent="parent"),
                CompoundInterval([0], [5], Strand.PLUS, parent="parent"),
                False,
                True,
            ),
            # Same coordinates different parent
            (
                CompoundInterval([0], [5], Strand.PLUS, parent="parent1"),
                CompoundInterval([0], [5], Strand.PLUS, parent="parent2"),
                False,
                False,
            ),
        ],
    )
    def test_contains_compound_interval(self, location1, location2, match_strand, expected):
        assert location1.contains(location2, match_strand) is expected

    @pytest.mark.parametrize(
        "single_interval,compound_interval,expected",
        [
            # Other location is exact interval, plus/minus strands
            (
                SingleInterval(3, 5, Strand.PLUS),
                CompoundInterval([3], [5], Strand.MINUS),
                SingleInterval(0, 2, Strand.MINUS),
            ),
            # One block of other location overlaps interval to the right, both have parents, minus/unstranded
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
            # One block of other location overlaps interval to the left, minus/plus
            (
                SingleInterval(10, 20, Strand.MINUS),
                CompoundInterval([5, 30], [12, 40], Strand.PLUS),
                SingleInterval(8, 10, Strand.MINUS),
            ),
            # Two blocks of other location overlap interval to the right, plus/minus
            (
                SingleInterval(3, 30, Strand.PLUS),
                CompoundInterval([10, 20, 40], [12, 35, 42], Strand.MINUS),
                CompoundInterval([7, 17], [9, 27], Strand.MINUS),
            ),
            # Two blocks of other location overlap interval to the left, plus/unstranded
            (
                SingleInterval(30, 50, Strand.PLUS),
                CompoundInterval([10, 25, 45], [20, 32, 47], Strand.UNSTRANDED),
                CompoundInterval([0, 15], [2, 17], Strand.UNSTRANDED),
            ),
            # Interval completely contains other location, minus/plus
            (
                SingleInterval(10, 30, Strand.MINUS),
                CompoundInterval([12, 22], [14, 24], Strand.PLUS),
                CompoundInterval([6, 16], [8, 18], Strand.MINUS),
            ),
            # Interval is completely contained in one block of other location, plus/plus
            (
                SingleInterval(3, 5, Strand.PLUS),
                CompoundInterval([0, 10], [8, 20], Strand.PLUS),
                SingleInterval(0, 2, Strand.PLUS),
            ),
        ],
    )
    def test_location_relative_to_single_interval(self, single_interval, compound_interval, expected):
        assert compound_interval.location_relative_to(single_interval) == expected

    @pytest.mark.parametrize(
        "single_interval,compound_interval,expected_exception",
        [
            # Interval is empty
            (
                SingleInterval(5, 5, Strand.PLUS),
                CompoundInterval([0], [10], Strand.PLUS),
                ValueError,
            ),
            # Other location is empty
            (
                SingleInterval(0, 10, Strand.PLUS),
                CompoundInterval([5], [5], Strand.PLUS),
                ValueError,
            ),
            # Interval has no parent, other location has a parent
            (
                SingleInterval(1, 2, Strand.PLUS),
                CompoundInterval([0], [5], Strand.PLUS, parent="parent"),
                ValueError,
            ),
            # Interval has a parent, other location doesn't
            (
                SingleInterval(1, 2, Strand.PLUS, parent="parent"),
                CompoundInterval([0], [5], Strand.PLUS),
                ValueError,
            ),
            # Parents are different
            (
                SingleInterval(1, 2, Strand.PLUS, parent="parent1"),
                CompoundInterval([0], [5], Strand.PLUS, parent="parent2"),
                ValueError,
            ),
            # No overlap
            (
                SingleInterval(0, 3, Strand.PLUS),
                CompoundInterval([5], [10], Strand.PLUS),
                ValueError,
            ),
        ],
    )
    def test_location_relative_to_error_single_interval(self, single_interval, compound_interval, expected_exception):
        with pytest.raises(expected_exception):
            compound_interval.location_relative_to(single_interval)

    @pytest.mark.parametrize(
        "location1,location2,expected",
        [
            # Both have parent
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS, parent="parent"),
                CompoundInterval([12, 20], [18, 25], Strand.MINUS, parent="parent"),
                SingleInterval(7, 10, Strand.MINUS),
            ),
            # One location is unstranded; overlaps two blocks
            (
                CompoundInterval([0, 10], [7, 17], Strand.PLUS),
                CompoundInterval([5, 20], [12, 30], Strand.UNSTRANDED),
                SingleInterval(5, 9, Strand.UNSTRANDED),
            ),
            # Location contains entire location; same strand
            (
                CompoundInterval([5, 10, 15], [7, 13, 19], Strand.MINUS),
                CompoundInterval([3, 8, 14], [8, 13, 20], Strand.MINUS),
                SingleInterval(0, 9, Strand.PLUS),
            ),
            # Locations equal
            (
                CompoundInterval([0, 3], [1, 10], Strand.PLUS),
                CompoundInterval([0, 3], [1, 10], Strand.PLUS),
                SingleInterval(0, 8, Strand.PLUS),
            ),
            # One location has multiple blocks within multiple blocks of other location
            (
                CompoundInterval([0, 20, 40], [10, 30, 50], Strand.MINUS),
                CompoundInterval([3, 24, 41], [22, 27, 45], Strand.PLUS),
                CompoundInterval([5, 13, 18], [9, 16, 27], Strand.MINUS),
            ),
        ],
    )
    def test_location_relative_to_compound_interval(self, location1, location2, expected):
        assert location2.location_relative_to(location1) == expected

    def test_location_relative_to_empty_location(self):
        assert CompoundInterval([5], [10], Strand.PLUS).location_relative_to(EmptyLocation()) is EmptyLocation()

    @pytest.mark.parametrize(
        "location1,location2,expected_exception",
        [
            # Different parents
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS, parent="parent1"),
                CompoundInterval([12, 20], [18, 25], Strand.MINUS, parent="parent2"),
                ValueError,
            ),
        ],
    )
    def test_location_relative_to_error_compound_interval(self, location1, location2, expected_exception):
        with pytest.raises(expected_exception):
            location2.location_relative_to(location1)

    @pytest.mark.parametrize(
        "location,expected_output",
        [
            (
                CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                CompoundLocation(
                    [
                        FeatureLocation(ExactPosition(0), ExactPosition(5), strand=1),
                        FeatureLocation(ExactPosition(10), ExactPosition(15), strand=1),
                    ],
                    "join",
                ),
            ),
            (
                CompoundInterval([0, 10], [5, 15], Strand.MINUS),
                CompoundLocation(
                    [
                        FeatureLocation(ExactPosition(0), ExactPosition(5), strand=-1),
                        FeatureLocation(ExactPosition(10), ExactPosition(15), strand=-1),
                    ],
                    "join",
                ),
            ),
            (
                CompoundInterval([0], [10], Strand.MINUS),
                FeatureLocation(ExactPosition(0), ExactPosition(10), strand=-1),
            ),
        ],
    )
    def test_biopython(self, location, expected_output):
        assert location.to_compound_location() == expected_output
