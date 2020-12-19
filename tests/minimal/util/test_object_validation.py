import pytest

from inscripta.biocantor.exc import (
    LocationOverlapException,
    LocationException,
    NullParentException,
    MismatchedParentException,
    NullSequenceException,
)
from inscripta.biocantor.location.location_impl import SingleInterval, CompoundInterval
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.parent import Parent
from inscripta.biocantor.sequence import Sequence
from inscripta.biocantor.util.object_validation import ObjectValidation


class TestObjectValidation:
    def test_require_location_nonempty(self):
        with pytest.raises(LocationException):
            ObjectValidation.require_location_nonempty(SingleInterval(5, 5, Strand.PLUS))
        with pytest.raises(LocationException):
            ObjectValidation.require_location_nonempty(CompoundInterval([5, 10], [5, 10], Strand.PLUS))
        ObjectValidation.require_location_nonempty(SingleInterval(5, 6, Strand.PLUS))

    def test_require_location_has_parent(self):
        with pytest.raises(NullParentException):
            ObjectValidation.require_location_has_parent(SingleInterval(0, 5, Strand.PLUS))
        ObjectValidation.require_location_has_parent(SingleInterval(0, 5, Strand.PLUS, parent="parent"))

    def test_require_location_has_parent_with_sequence(self):
        with pytest.raises(NullParentException):
            ObjectValidation.require_location_has_parent_with_sequence(SingleInterval(0, 5, Strand.PLUS))
        with pytest.raises(NullSequenceException):
            ObjectValidation.require_location_has_parent_with_sequence(
                SingleInterval(0, 5, Strand.PLUS, parent="parent")
            )
        ObjectValidation.require_location_has_parent_with_sequence(
            SingleInterval(
                0,
                5,
                Strand.PLUS,
                parent=Parent(id="parent", sequence=Sequence("AAAAA", Alphabet.NT_STRICT)),
            )
        )

    def test_require_parent_has_location(self):
        with pytest.raises(NullParentException):
            ObjectValidation.require_parent_has_location(Parent(id="parent"))
        ObjectValidation.require_parent_has_location(Parent(location=SingleInterval(5, 6, Strand.PLUS)))

    def test_require_parent_has_parent(self):
        with pytest.raises(NullParentException):
            ObjectValidation.require_parent_has_parent(Parent(id="parent"))
        ObjectValidation.require_parent_has_parent(Parent(id="parent", parent="grandparent"))

    def test_require_parent_has_parent_with_location(self):
        with pytest.raises(NullParentException):
            ObjectValidation.require_parent_has_parent_with_location(Parent(id="parent"))
        with pytest.raises(NullParentException):
            ObjectValidation.require_parent_has_parent_with_location(Parent(id="parent", parent="grandparent"))
        ObjectValidation.require_parent_has_parent_with_location(
            Parent(
                id="parent",
                parent=Parent(id="grandparent", location=SingleInterval(0, 5, Strand.PLUS)),
            )
        )

    def test_require_parents_equal_except_location(self):
        with pytest.raises(MismatchedParentException):
            ObjectValidation.require_parents_equal_except_location(Parent(id="parent1"), Parent(id="parent2"))
        ObjectValidation.require_parents_equal_except_location(
            Parent(id="parent", location=SingleInterval(0, 5, Strand.PLUS)),
            Parent(id="parent", location=SingleInterval(10, 20, Strand.MINUS)),
        )

    def require_locations_have_same_nonempty_parent(self):
        raise NotImplementedError

    def test_require_locations_overlap(self):
        with pytest.raises(LocationOverlapException):
            ObjectValidation.require_locations_overlap(
                SingleInterval(0, 5, Strand.PLUS), SingleInterval(5, 10, Strand.PLUS)
            )
        with pytest.raises(LocationOverlapException):
            ObjectValidation.require_locations_overlap(
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(0, 5, Strand.MINUS),
                match_strand=True,
            )
        ObjectValidation.require_locations_overlap(
            SingleInterval(0, 5, Strand.PLUS),
            SingleInterval(3, 6, Strand.PLUS),
            match_strand=True,
        )

    def test_require_locations_do_not_overlap(self):
        with pytest.raises(LocationOverlapException):
            ObjectValidation.require_locations_do_not_overlap(
                SingleInterval(0, 5, Strand.PLUS), SingleInterval(3, 8, Strand.PLUS)
            )
        with pytest.raises(LocationOverlapException):
            ObjectValidation.require_locations_do_not_overlap(
                SingleInterval(0, 5, Strand.PLUS),
                SingleInterval(3, 8, Strand.MINUS),
                match_strand=False,
            )
        ObjectValidation.require_locations_do_not_overlap(
            SingleInterval(0, 5, Strand.PLUS),
            SingleInterval(5, 10, Strand.PLUS),
            match_strand=True,
        )

    def test_require_object_has_type(self):
        with pytest.raises(TypeError):
            ObjectValidation.require_object_has_type(
                CompoundInterval.from_single_intervals([SingleInterval(5, 10, Strand.PLUS)]),
                SingleInterval,
            )
        ObjectValidation.require_object_has_type(
            CompoundInterval.from_single_intervals([SingleInterval(5, 10, Strand.PLUS)]),
            CompoundInterval,
        )
