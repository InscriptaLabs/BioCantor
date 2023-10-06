import pytest

from biocantor.exc import (
    NoSuchAncestorException,
    InvalidStrandException,
    InvalidPositionException,
    NullParentException,
    MismatchedParentException,
    ParentException,
    LocationException,
)
from biocantor.location.location_impl import (
    SingleInterval,
    CompoundInterval,
    EmptyLocation,
)
from biocantor.location.strand import Strand
from biocantor.parent import Parent, make_parent
from biocantor.sequence.alphabet import Alphabet
from biocantor.sequence import Sequence


class TestParent:
    @pytest.mark.parametrize(
        "id,sequence_type,strand,location,sequence,expected",
        [
            (None, None, None, None, None, Parent()),
            (
                "id",
                "seqtype",
                Strand.MINUS,
                SingleInterval(
                    0,
                    1,
                    Strand.MINUS,
                    Parent(
                        id="id",
                        sequence_type="seqtype",
                        strand=Strand.MINUS,
                        sequence=Sequence(
                            "AAA",
                            Alphabet.NT_STRICT,
                            id="id",
                            type="seqtype",
                            parent=Parent(
                                id="id2",
                                sequence=Sequence("CCC", Alphabet.NT_STRICT, id="id2"),
                            ),
                        ),
                    ),
                ),
                Sequence(
                    "AAA",
                    Alphabet.NT_STRICT,
                    id="id",
                    type="seqtype",
                    parent=Parent(id="id2", sequence=Sequence("CCC", Alphabet.NT_STRICT, id="id2")),
                ),
                Parent(
                    id="id",
                    sequence_type="seqtype",
                    strand=Strand.MINUS,
                    location=SingleInterval(
                        0,
                        1,
                        Strand.MINUS,
                        Parent(
                            id="id",
                            sequence_type="seqtype",
                            strand=Strand.MINUS,
                            sequence=Sequence(
                                "AAA",
                                Alphabet.NT_STRICT,
                                id="id",
                                type="seqtype",
                                parent=Parent(
                                    id="id2",
                                    sequence=Sequence("CCC", Alphabet.NT_STRICT, id="id2"),
                                ),
                            ),
                        ),
                    ),
                    sequence=Sequence(
                        "AAA",
                        Alphabet.NT_STRICT,
                        id="id",
                        type="seqtype",
                        parent=Parent(
                            id="id2",
                            sequence=Sequence("CCC", Alphabet.NT_STRICT, id="id2"),
                        ),
                    ),
                ),
            ),
        ],
    )
    def test_init(self, id, sequence_type, strand, location, sequence, expected):
        assert expected == Parent(
            id=id,
            sequence_type=sequence_type,
            strand=strand,
            location=location,
            sequence=sequence,
        )

    @pytest.mark.parametrize(
        "id,sequence_type,strand,location,sequence,parent,expected_exception",
        [
            ("id1", None, None, SingleInterval(0, 5, Strand.PLUS, parent="id2"), None, None, ParentException),
            ("id1", None, None, None, Sequence("AAA", Alphabet.NT_STRICT, id="id2"), None, ParentException),
            (
                None,
                None,
                None,
                SingleInterval(0, 5, Strand.PLUS, parent="id1"),
                Sequence("AAC", Alphabet.NT_STRICT, id="id2"),
                None,
                ParentException,
            ),
            (
                None,
                "seqtype",
                None,
                SingleInterval(0, 5, Strand.PLUS, parent=Parent(sequence_type="unknown")),
                None,
                None,
                ParentException,
            ),
            (None, "seqtype", None, None, Sequence("AAT", Alphabet.NT_STRICT, type="unknown"), None, ParentException),
            (
                None,
                None,
                None,
                SingleInterval(0, 5, Strand.PLUS, parent=Parent(sequence_type="unknown")),
                Sequence("AAG", Alphabet.NT_STRICT, type="seqtype"),
                None,
                ParentException,
            ),
            (None, None, Strand.MINUS, SingleInterval(0, 5, Strand.PLUS), None, None, InvalidStrandException),
            (
                None,
                None,
                None,
                SingleInterval(0, 10, Strand.PLUS),
                Sequence("A", Alphabet.NT_STRICT),
                None,
                InvalidPositionException,
            ),
            (
                None,
                None,
                None,
                None,
                Sequence("AA", Alphabet.NT_STRICT),
                Parent(sequence=Sequence("A", Alphabet.NT_STRICT)),
                LocationException,
            ),
            (None, None, Strand.PLUS, SingleInterval(5, 10, Strand.MINUS), None, None, InvalidStrandException),
            (
                None,
                None,
                None,
                None,
                Sequence("AA", Alphabet.NT_STRICT, parent="id1"),
                Parent(id="id2"),
                MismatchedParentException,
            ),
        ],
    )
    def test_init_error(self, id, sequence_type, strand, location, sequence, parent, expected_exception):
        with pytest.raises(expected_exception):
            Parent(
                id=id,
                sequence_type=sequence_type,
                strand=strand,
                location=location,
                sequence=sequence,
                parent=parent,
            )

    @pytest.mark.parametrize(
        "obj,expected",
        [
            (
                Sequence("AAA", Alphabet.NT_STRICT),
                Parent(sequence=Sequence("AAA", Alphabet.NT_STRICT)),
            ),
            ("parent", Parent(id="parent")),
            (
                SingleInterval(5, 10, Strand.PLUS),
                Parent(location=SingleInterval(5, 10, Strand.PLUS)),
            ),
            (
                CompoundInterval([5], [10], Strand.PLUS),
                Parent(location=CompoundInterval([5], [10], Strand.PLUS)),
            ),
            (EmptyLocation(), Parent(location=EmptyLocation())),
            (Strand.MINUS, Parent(strand=Strand.MINUS)),
            (
                Parent(
                    id="parent",
                    sequence_type="chr",
                    strand=Strand.MINUS,
                    parent=Parent(id="grandparent"),
                ),
                Parent(
                    id="parent",
                    sequence_type="chr",
                    strand=Strand.MINUS,
                    parent=Parent(id="grandparent"),
                ),
            ),
        ],
    )
    def test_make_parent(self, obj, expected):
        assert make_parent(obj) == expected

    @pytest.mark.parametrize(
        "parent1,parent2,expected",
        [
            (Parent(), Parent(), True),
            (Parent(), Parent(id=None, sequence_type=None), True),
            (Parent(id="id1"), Parent(id="id2"), False),
            (
                Parent(sequence_type=None),
                Parent(sequence_type="unknown"),
                False,
            ),
            (Parent(strand=Strand.UNSTRANDED), Parent(strand=Strand.MINUS), False),
            (
                Parent(location=SingleInterval(0, 5, Strand.PLUS, parent="id1")),
                Parent(location=SingleInterval(0, 5, Strand.PLUS, parent="id2")),
                False,
            ),
            (
                Parent(sequence=Sequence("A", Alphabet.NT_STRICT)),
                Parent(sequence=Sequence("A", Alphabet.NT_STRICT, parent=Parent(id="parent"))),
                False,
            ),
            (
                Parent(parent="parent1"),
                Parent(parent="parent2"),
                False,
            ),
        ],
    )
    def test_eq(self, parent1, parent2, expected):
        assert (parent1 == parent2) is expected

    @pytest.mark.parametrize(
        "parent1,parent2,expected",
        [
            (Parent(), Parent(), True),
            (Parent(id="id1"), Parent(id="id2"), False),
            (
                Parent(sequence_type=None),
                Parent(sequence_type="unknown"),
                False,
            ),
            (Parent(strand=Strand.UNSTRANDED), Parent(strand=Strand.MINUS), True),
            (
                Parent(location=SingleInterval(0, 5, Strand.PLUS, parent="id1")),
                Parent(location=SingleInterval(0, 5, Strand.PLUS, parent="id2")),
                False,
            ),
            (
                Parent(sequence=Sequence("A", Alphabet.NT_STRICT)),
                Parent(sequence=Sequence("A", Alphabet.NT_STRICT, parent="parent")),
                False,
            ),
            (
                Parent(parent="parent1"),
                Parent(parent="parent2"),
                False,
            ),
        ],
    )
    def test_equals_except_location(self, parent1, parent2, expected):
        assert parent1.equals_except_location(parent2) is expected

    @pytest.mark.parametrize(
        "id,location,sequence,expected",
        [
            ("id", None, None, "id"),
            (
                None,
                SingleInterval(0, 1, Strand.PLUS, parent="id"),
                None,
                "id",
            ),
            (
                None,
                None,
                Sequence("A", Alphabet.NT_STRICT, id="id", parent="id2"),
                "id",
            ),
            (
                "id",
                SingleInterval(0, 1, Strand.PLUS, parent="id"),
                Sequence("A", Alphabet.NT_STRICT, id="id", parent="id2"),
                "id",
            ),
        ],
    )
    def test_id(self, id, location, sequence, expected):
        assert Parent(id=id, location=location, sequence=sequence).id == expected

    @pytest.mark.parametrize(
        "sequence_type,location,sequence,expected",
        [
            ("seqtype", None, None, "seqtype"),
            (
                None,
                SingleInterval(
                    0,
                    5,
                    Strand.PLUS,
                    parent=Parent(sequence_type="seqtype"),
                ),
                None,
                "seqtype",
            ),
            (
                None,
                None,
                Sequence("A", Alphabet.NT_STRICT, type="seqtype"),
                "seqtype",
            ),
            (
                None,
                None,
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    type="seqtype",
                    parent=Parent(sequence_type="seqtype_2"),
                ),
                "seqtype",
            ),
        ],
    )
    def test_sequence_type(self, sequence_type, location, sequence, expected):
        assert Parent(sequence_type=sequence_type, location=location, sequence=sequence).sequence_type == expected

    @pytest.mark.parametrize(
        "strand,location,sequence,expected",
        [
            (Strand.PLUS, None, None, Strand.PLUS),
            (None, SingleInterval(0, 5, Strand.MINUS), None, Strand.MINUS),
            (
                Strand.PLUS,
                None,
                Sequence("A", Alphabet.NT_STRICT, parent=Strand.MINUS),
                Strand.PLUS,
            ),
        ],
    )
    def test_strand(self, strand, location, sequence, expected):
        assert Parent(strand=strand, location=location, sequence=sequence).strand == expected

    def test_location(self):
        assert Parent(location=SingleInterval(0, 1, Strand.PLUS)).location == SingleInterval(0, 1, Strand.PLUS)

    def test_sequence(self):
        assert Parent(sequence=Sequence("A", Alphabet.NT_STRICT)).sequence == Sequence("A", Alphabet.NT_STRICT)

    @pytest.mark.parametrize(
        "parent,expected",
        [
            (Parent(parent="id"), Parent(id="id")),
            (
                Parent(
                    sequence=Sequence(
                        "AA",
                        Alphabet.NT_STRICT,
                        parent=Parent(sequence_type="chr"),
                    )
                ),
                Parent(sequence_type="chr"),
            ),
        ],
    )
    def test_parent(self, parent, expected):
        assert parent.parent == expected

    @pytest.mark.parametrize(
        "parent,expected",
        [
            (Parent(), Parent()),
            (Parent(strand=Strand.PLUS), Parent()),
            (
                Parent(strand=Strand.PLUS, location=SingleInterval(5, 10, Strand.PLUS)),
                Parent(),
            ),
            (
                Parent(
                    id="parent",
                    sequence_type="unknown",
                    strand=Strand.PLUS,
                    location=SingleInterval(0, 1, Strand.PLUS),
                    sequence=Sequence("AAA", Alphabet.NT_STRICT),
                    parent="grandparent",
                ),
                Parent(
                    id="parent",
                    sequence_type="unknown",
                    sequence=Sequence("AAA", Alphabet.NT_STRICT),
                    parent="grandparent",
                ),
            ),
        ],
    )
    def test_strip_location_info(self, parent, expected):
        assert parent.strip_location_info() == expected

    @pytest.mark.parametrize(
        "parent,sequence_type,include_self,expected",
        [
            (
                Parent(
                    id="self",
                    sequence_type="seqtype",
                    parent=Parent(id="parent", sequence_type="seqtype"),
                ),
                "seqtype",
                True,
                Parent(
                    id="self",
                    sequence_type="seqtype",
                    parent=Parent(id="parent", sequence_type="seqtype"),
                ),
            ),
            (
                Parent(
                    id="self",
                    sequence_type="seqtype",
                    parent=Parent(id="parent", sequence_type="seqtype"),
                ),
                "seqtype",
                False,
                Parent(id="parent", sequence_type="seqtype"),
            ),
            (
                Parent(
                    id="self",
                    sequence_type="seqtype",
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype_2",
                        parent=Parent(id="grandparent", sequence_type="seqtype_2"),
                    ),
                ),
                "seqtype_2",
                True,
                Parent(
                    id="parent",
                    sequence_type="seqtype_2",
                    parent=Parent(id="grandparent", sequence_type="seqtype_2"),
                ),
            ),
        ],
    )
    def test_first_ancestor_of_type(self, parent, sequence_type, include_self, expected):
        assert parent.first_ancestor_of_type(sequence_type, include_self=include_self) == expected

    @pytest.mark.parametrize(
        "parent,sequence_type,include_self",
        [
            (Parent(id="self"), "seqtype_2", True),
            (
                Parent(id="self", parent="parent"),
                "seqtype_2",
                True,
            ),
            (
                Parent(
                    id="self",
                    sequence_type="seqtype",
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype_2",
                        parent=Parent(id="grandparent", sequence_type="seqtype_2"),
                    ),
                ),
                "chr",
                True,
            ),
        ],
    )
    def test_first_ancestor_of_type_error(self, parent, sequence_type, include_self):
        with pytest.raises(NoSuchAncestorException):
            parent.first_ancestor_of_type(sequence_type, include_self=include_self)

    @pytest.mark.parametrize(
        "parent,sequence_type,include_self,expected",
        [
            (
                Parent(
                    id="self",
                    sequence_type="seqtype",
                    parent=Parent(id="parent", sequence_type="seqtype"),
                ),
                "seqtype",
                True,
                True,
            ),
            (
                Parent(
                    id="self",
                    sequence_type="seqtype",
                    parent=Parent(id="parent", sequence_type="seqtype"),
                ),
                "seqtype",
                False,
                True,
            ),
            (
                Parent(
                    id="self",
                    sequence_type="seqtype",
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype_2",
                        parent=Parent(id="grandparent", sequence_type="seqtype_2"),
                    ),
                ),
                "seqtype_2",
                True,
                True,
            ),
            (
                Parent(id="self"),
                "seqtype_2",
                True,
                False,
            ),
            (
                Parent(id="self", parent="parent"),
                "seqtype_2",
                True,
                False,
            ),
            (
                Parent(
                    id="self",
                    sequence_type="seqtype",
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype_2",
                        parent=Parent(id="grandparent", sequence_type="seqtype_2"),
                    ),
                ),
                "chr",
                True,
                False,
            ),
        ],
    )
    def test_has_ancestor_of_type(self, parent, sequence_type, include_self, expected):
        assert parent.has_ancestor_of_type(sequence_type, include_self=include_self) is expected

    @pytest.mark.parametrize(
        "parent,expected",
        [
            (
                Parent(
                    id="parent",
                    location=SingleInterval(3, 5, Strand.PLUS),
                    parent=Parent(id="grandparent", location=SingleInterval(10, 20, Strand.PLUS)),
                ),
                SingleInterval(13, 15, Strand.PLUS, parent="grandparent"),
            ),
            (
                Parent(
                    id="parent",
                    location=SingleInterval(0, 5, Strand.PLUS),
                    sequence_type="unknown",
                    strand=Strand.PLUS,
                    parent=Parent(
                        id="grandparent",
                        location=SingleInterval(100, 200, Strand.MINUS),
                    ),
                ),
                SingleInterval(195, 200, Strand.MINUS, parent="grandparent"),
            ),
            (
                Parent(
                    id="parent",
                    location=SingleInterval(6, 9, Strand.MINUS),
                    parent=Parent(
                        id="grandparent",
                        location=SingleInterval(0, 10, Strand.PLUS),
                        sequence_type="chr",
                        strand=Strand.PLUS,
                        parent="great grandparent",
                    ),
                ),
                SingleInterval(
                    6,
                    9,
                    Strand.MINUS,
                    parent=Parent(
                        id="grandparent",
                        sequence_type="chr",
                        parent="great grandparent",
                    ),
                ),
            ),
            (
                Parent(
                    id="parent",
                    sequence_type="chr",
                    strand=Strand.MINUS,
                    location=SingleInterval(6, 8, Strand.MINUS),
                    parent=Parent(
                        id="grandparent",
                        sequence_type="unknown",
                        strand=Strand.MINUS,
                        location=SingleInterval(5, 15, Strand.MINUS),
                        parent="great grandparent",
                    ),
                ),
                SingleInterval(
                    7,
                    9,
                    Strand.PLUS,
                    parent=Parent(
                        id="grandparent",
                        sequence_type="unknown",
                        parent="great grandparent",
                    ),
                ),
            ),
            (
                Parent(
                    id="parent",
                    location=SingleInterval(3, 5, Strand.UNSTRANDED),
                    parent=Parent(id="grandparent", location=SingleInterval(10, 20, Strand.PLUS)),
                ),
                SingleInterval(13, 15, Strand.UNSTRANDED, parent="grandparent"),
            ),
            (
                Parent(
                    id="parent",
                    location=SingleInterval(3, 5, Strand.UNSTRANDED),
                    parent=Parent(id="grandparent", location=SingleInterval(10, 20, Strand.MINUS)),
                ),
                SingleInterval(15, 17, Strand.UNSTRANDED, parent="grandparent"),
            ),
        ],
    )
    def test_lift_child_location_contiguous_to_parent_single_interval(self, parent, expected):
        assert parent.lift_child_location_to_parent() == expected

    @pytest.mark.parametrize(
        "parent,expected",
        [
            (
                Parent(
                    id="parent",
                    location=CompoundInterval([3, 7], [5, 10], Strand.PLUS),
                    parent=Parent(id="grandparent", location=SingleInterval(10, 20, Strand.PLUS)),
                ),
                CompoundInterval([13, 17], [15, 20], Strand.PLUS, parent="grandparent"),
            ),
            (
                Parent(
                    id="parent",
                    location=CompoundInterval([0, 10], [5, 15], Strand.PLUS),
                    sequence_type="unknown",
                    strand=Strand.PLUS,
                    parent=Parent(
                        id="grandparent",
                        location=SingleInterval(100, 200, Strand.MINUS),
                    ),
                ),
                CompoundInterval(
                    [185, 195],
                    [190, 200],
                    Strand.MINUS,
                    parent="grandparent",
                ),
            ),
            (
                Parent(
                    id="parent",
                    location=CompoundInterval([6], [9], Strand.MINUS),
                    parent=Parent(
                        id="grandparent",
                        location=SingleInterval(0, 10, Strand.PLUS),
                        sequence_type="chr",
                        strand=Strand.PLUS,
                        parent="great grandparent",
                    ),
                ),
                SingleInterval(
                    6,
                    9,
                    Strand.MINUS,
                    parent=Parent(
                        id="grandparent",
                        sequence_type="chr",
                        parent="great grandparent",
                    ),
                ),
            ),
            (
                Parent(
                    id="parent",
                    sequence_type="chr",
                    strand=Strand.MINUS,
                    location=CompoundInterval([6], [8], Strand.MINUS),
                    parent=Parent(
                        id="grandparent",
                        sequence_type="unknown",
                        strand=Strand.MINUS,
                        location=SingleInterval(5, 15, Strand.MINUS),
                        parent="great grandparent",
                    ),
                ),
                SingleInterval(
                    7,
                    9,
                    Strand.PLUS,
                    parent=Parent(
                        id="grandparent",
                        sequence_type="unknown",
                        parent="great grandparent",
                    ),
                ),
            ),
            (
                Parent(
                    id="parent",
                    location=CompoundInterval([3, 7], [5, 10], Strand.UNSTRANDED),
                    parent=Parent(id="grandparent", location=SingleInterval(10, 20, Strand.PLUS)),
                ),
                CompoundInterval(
                    [13, 17],
                    [15, 20],
                    Strand.UNSTRANDED,
                    parent="grandparent",
                ),
            ),
            (
                Parent(
                    id="parent",
                    location=CompoundInterval([3], [5], Strand.UNSTRANDED),
                    parent=Parent(id="grandparent", location=SingleInterval(10, 20, Strand.MINUS)),
                ),
                SingleInterval(15, 17, Strand.UNSTRANDED, parent="grandparent"),
            ),
        ],
    )
    def test_lift_child_location_discontiguous_to_parent_single_interval(self, parent, expected):
        assert parent.lift_child_location_to_parent() == expected

    @pytest.mark.parametrize(
        "parent,expected_error",
        [
            # No location
            (
                Parent(parent=SingleInterval(5, 10, Strand.PLUS)),
                NullParentException,
            ),
            # Parent has no location
            (
                Parent(
                    location=SingleInterval(5, 10, Strand.PLUS),
                    parent="grandparent",
                ),
                NullParentException,
            ),
            # Location on parent can't be unstranded
            (
                Parent(
                    location=SingleInterval(5, 10, Strand.PLUS),
                    parent=Parent(
                        id="grandparent",
                        location=SingleInterval(0, 100, Strand.UNSTRANDED),
                    ),
                ),
                InvalidStrandException,
            ),
            # Location must fit inside location on parent
            (
                Parent(
                    location=SingleInterval(5, 10, Strand.PLUS),
                    parent=Parent(id="grandparent", location=SingleInterval(30, 31, Strand.PLUS)),
                ),
                ValueError,
            ),
        ],
    )
    def test_lift_child_location_to_parent_single_interval_error(self, parent, expected_error):
        with pytest.raises(expected_error):
            parent.lift_child_location_to_parent()

    @pytest.mark.parametrize(
        "parent,expected",
        [
            # Location takes up entire parent location
            (
                Parent(
                    id="parent",
                    location=SingleInterval(0, 10, Strand.PLUS),
                    parent=Parent(
                        id="grandparent",
                        location=CompoundInterval([5, 20], [10, 25], Strand.PLUS),
                    ),
                ),
                CompoundInterval([5, 20], [10, 25], Strand.PLUS, parent="grandparent"),
            ),
            # Location (unstranded) takes up part of parent location (minus)
            (
                Parent(
                    id="parent",
                    location=SingleInterval(10, 20, Strand.UNSTRANDED),
                    parent=Parent(
                        id="grandparent",
                        location=CompoundInterval([10, 20, 30], [18, 28, 38], Strand.MINUS),
                    ),
                ),
                CompoundInterval(
                    [14, 20],
                    [18, 26],
                    Strand.UNSTRANDED,
                    parent="grandparent",
                ),
            ),
            # Location (minus) takes up one block of parent location (plus); location is at end of sequence
            (
                Parent(
                    id="parent",
                    location=SingleInterval(5, 10, Strand.MINUS),
                    parent=Parent(
                        id="grandparent",
                        location=CompoundInterval([30, 40], [35, 45], Strand.PLUS),
                    ),
                ),
                SingleInterval(40, 45, Strand.MINUS, parent="grandparent"),
            ),
            # Location (minus) takes up part of one block of parent location (minus)
            (
                Parent(
                    id="parent",
                    location=SingleInterval(0, 4, Strand.MINUS),
                    parent=Parent(
                        id="grandparent",
                        location=CompoundInterval([30, 40], [35, 45], Strand.MINUS),
                    ),
                ),
                SingleInterval(41, 45, Strand.PLUS, parent="grandparent"),
            ),
        ],
    )
    def test_lift_child_location_contiguous_to_parent_compound_interval(self, parent, expected):
        assert parent.lift_child_location_to_parent() == expected

    @pytest.mark.parametrize(
        "parent,expected",
        [
            # Location takes up entire parent location
            (
                Parent(
                    id="parent",
                    location=CompoundInterval([0, 5], [5, 10], Strand.PLUS),
                    parent=Parent(
                        id="grandparent",
                        location=CompoundInterval([5, 20], [10, 25], Strand.PLUS),
                    ),
                ),
                CompoundInterval([5, 20], [10, 25], Strand.PLUS, parent="grandparent"),
            ),
            # Location (unstranded) takes up part of parent location (minus)
            (
                Parent(
                    id="parent",
                    location=CompoundInterval([10, 22], [20, 23], Strand.UNSTRANDED),
                    parent=Parent(
                        id="grandparent",
                        location=CompoundInterval([10, 20, 30], [18, 28, 38], Strand.MINUS),
                    ),
                ),
                CompoundInterval(
                    [11, 14, 20],
                    [12, 18, 26],
                    Strand.UNSTRANDED,
                    parent="grandparent",
                ),
            ),
            # Location (minus) takes up one block of parent location (plus); location is at end of sequence
            (
                Parent(
                    id="parent",
                    location=CompoundInterval([5], [10], Strand.MINUS),
                    parent=Parent(
                        id="grandparent",
                        location=CompoundInterval([30, 40], [35, 45], Strand.PLUS),
                    ),
                ),
                SingleInterval(40, 45, Strand.MINUS, parent="grandparent"),
            ),
            # Location (minus) takes up part of one block of parent location (minus)
            (
                Parent(
                    id="parent",
                    location=CompoundInterval([0, 3], [1, 4], Strand.MINUS),
                    parent=Parent(
                        id="grandparent",
                        location=CompoundInterval([30, 40], [35, 45], Strand.MINUS),
                    ),
                ),
                CompoundInterval([41, 44], [42, 45], Strand.PLUS, parent="grandparent"),
            ),
        ],
    )
    def test_lift_child_location_discontiguous_to_parent_compound_interval(self, parent, expected):
        assert parent.lift_child_location_to_parent() == expected

    @pytest.mark.parametrize(
        "parent,expected_error",
        [
            # Location must fit inside location on parent
            (
                Parent(
                    location=SingleInterval(5, 50, Strand.PLUS),
                    parent=Parent(
                        id="grandparent",
                        location=CompoundInterval([10, 20], [15, 25], Strand.PLUS),
                    ),
                ),
                InvalidPositionException,
            ),
        ],
    )
    def test_lift_child_location_to_parent_compound_interval_error(self, parent, expected_error):
        with pytest.raises(expected_error):
            parent.lift_child_location_to_parent()

    @pytest.mark.parametrize(
        "parent,location,expected",
        [
            (
                Parent(),
                SingleInterval(5, 10, Strand.PLUS),
                Parent(location=SingleInterval(5, 10, Strand.PLUS), strand=Strand.PLUS),
            ),
            (
                Parent(
                    id="parent",
                    sequence_type="unknown",
                    strand=Strand.MINUS,
                    location=SingleInterval(0, 2, Strand.MINUS),
                    sequence=Sequence("AAA", Alphabet.NT_STRICT),
                ),
                SingleInterval(2, 3, Strand.PLUS),
                Parent(
                    id="parent",
                    sequence_type="unknown",
                    strand=Strand.PLUS,
                    location=SingleInterval(2, 3, Strand.PLUS),
                    sequence=Sequence("AAA", Alphabet.NT_STRICT),
                ),
            ),
            (
                Parent(
                    id="parent",
                    sequence_type="unknown",
                    strand=Strand.MINUS,
                    location=SingleInterval(0, 2, Strand.MINUS),
                    sequence=Sequence("AAA", Alphabet.NT_STRICT),
                ),
                None,
                Parent(
                    id="parent",
                    sequence_type="unknown",
                    sequence=Sequence("AAA", Alphabet.NT_STRICT),
                ),
            ),
        ],
    )
    def test_reset_location(self, parent, location, expected):
        assert parent.reset_location(location) == expected

    @pytest.mark.parametrize(
        "parent,location,expected_exception",
        [
            (
                Parent(sequence=Sequence("AAA", Alphabet.NT_STRICT)),
                SingleInterval(0, 5, Strand.PLUS),
                InvalidPositionException,
            ),
            (
                Parent(id="id1", sequence=Sequence("AAA", Alphabet.NT_STRICT)),
                SingleInterval(
                    0,
                    1,
                    Strand.PLUS,
                    parent=Parent(id="id2", sequence=Sequence("AAA", Alphabet.NT_STRICT)),
                ),
                ParentException,
            ),
        ],
    )
    def test_reset_location_error(self, parent, location, expected_exception):
        with pytest.raises(expected_exception):
            parent.reset_location(location)

    @pytest.mark.parametrize(
        "parent,sequence,include_self,expected",
        [
            (Parent(), Sequence("AA", Alphabet.NT_STRICT), True, False),
            (Parent(), Sequence("AA", Alphabet.NT_STRICT), False, False),
            (
                Parent(sequence=Sequence("AA", Alphabet.NT_STRICT)),
                Sequence("AA", Alphabet.NT_STRICT),
                True,
                True,
            ),
            (
                Parent(sequence=Sequence("AA", Alphabet.NT_STRICT)),
                Sequence("AA", Alphabet.NT_STRICT),
                False,
                False,
            ),
            (
                Parent(
                    sequence=Sequence("AA", Alphabet.NT_STRICT),
                    parent=Sequence("AA", Alphabet.NT_STRICT),
                ),
                Sequence("AA", Alphabet.NT_STRICT),
                False,
                True,
            ),
            (
                Parent(
                    sequence=Sequence("AA", Alphabet.NT_STRICT),
                    parent=Sequence("AAT", Alphabet.NT_STRICT),
                ),
                Sequence("AAT", Alphabet.NT_STRICT),
                False,
                True,
            ),
            (
                Parent(
                    sequence=Sequence("AA", Alphabet.NT_STRICT),
                    parent=Sequence("AAT", Alphabet.NT_STRICT),
                ),
                Sequence("AAT", Alphabet.NT_STRICT, id="id"),
                True,
                False,
            ),
            (
                Parent(
                    parent=Parent(parent=Parent(parent=Parent(sequence=Sequence("AAA", Alphabet.NT_STRICT, id="seq"))))
                ),
                Sequence("AAA", Alphabet.NT_STRICT, id="seq"),
                True,
                True,
            ),
        ],
    )
    def test_has_ancestor_sequence(self, parent, sequence, include_self, expected):
        assert parent.has_ancestor_sequence(sequence, include_self) == expected
