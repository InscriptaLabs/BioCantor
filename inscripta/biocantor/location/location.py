from abc import ABC
from typing import Iterator, Union

from inscripta.biocantor import AbstractLocation
from inscripta.biocantor.exc import NoSuchAncestorException, NullParentException
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import Parent, SequenceType
from inscripta.biocantor.sequence import Sequence
from inscripta.biocantor.util.object_validation import ObjectValidation


class Location(AbstractLocation, ABC):
    """Abstract location with respect to a coordinate system"""

    def parent_to_relative_location(self, parent_location: "Location", optimize_blocks: bool = True) -> "Location":
        """Converts a Location on the parent to a Location relative to this Location.

        Parameters
        ----------
        parent_location
            Location with the same parent as this Location. Both parents can be None.
        optimize_blocks
            Run optimize_blocks on the resulting location?

        Returns
        -------
        New Location relative to this Location.
        """
        return parent_location.location_relative_to(self, optimize_blocks=optimize_blocks)

    def location_relative_to(self, other: "Location", optimize_blocks: bool = True) -> "Location":
        """Converts this Location to a Location relative to another Location. The Locations must overlap.
        The returned value represents the relative location of the overlap within the other Location.

        If ``optimize_blocks`` is ``True``, the resulting Location will not have any adjacent or overlapping
        intervals. This is often desirable, because the output of this function can have weird coordinates
        when the locations are overlapping or adjacent. However, there are some cases where it is desirable
        to retain the original block structure. One such example are CDS where adjacent blocks or overlapping
        blocks are used to model frameshifts or indels.
        """
        if other.is_empty:
            return other
        if self.parent or other.parent:
            if not self.parent and other.parent:
                raise NullParentException(
                    "Parents must be both null or both non-null:\n{}\n  !=\n{}".format(self.parent, other.parent)
                )
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        ObjectValidation.require_locations_overlap(self, other)
        return self._location_relative_to(other, optimize_blocks=optimize_blocks)

    def scan_windows(self, window_size: int, step_size: int, start_pos: int = 0) -> Iterator["Location"]:
        """Returns an iterator over fixed size windows within this Location. Windows represent sub-regions
        of this Location and are with respect to the same parent as this Location. The final window returned
        is the last one that fits completely within this Location. Returned windows are in order according to
        relative position within this Location; i.e., corresponding to the strand of this Location.

        Parameters
        ----------
        window_size
        step_size
        start_pos
            0-based relative start position of first window relative to this Location
        """
        if not 0 <= start_pos < len(self):
            raise ValueError("Start position ({}) must be within location length ({})".format(start_pos, len(self)))
        if min(window_size, step_size) < 1:
            raise ValueError("Window size and step size must both be positive: {}, {}".format(window_size, step_size))
        if window_size > len(self):
            raise ValueError("Window size ({}) must be <= location length ({})".format(window_size, len(self)))
        if start_pos + window_size > len(self):
            raise ValueError(
                "Start position ({}) + window size ({}) must be <= location length ({})".format(
                    start_pos, window_size, len(self)
                )
            )
        self.strand.assert_directional()
        for curr_start in range(start_pos, len(self) - window_size + 1, step_size):
            yield self.relative_interval_to_parent_location(curr_start, curr_start + window_size, Strand.PLUS)

    def first_ancestor_of_type(self, sequence_type: Union[str, SequenceType]) -> Parent:
        """Returns the Parent object representing the closest ancestor (parent, parent of parent, etc.)
        of this location which has the given sequence type. Raises NoSuchAncestorException if no ancestor with
        the given type exists."""
        if not self.parent:
            raise NoSuchAncestorException("Location has no parent")
        return self.parent.first_ancestor_of_type(sequence_type, include_self=True)

    def has_ancestor_of_type(self, sequence_type: Union[str, SequenceType]) -> bool:
        """Returns True if some ancestor (parent, parent of parent, etc.) of of this location has the given sequence
        type, or False otherwise."""
        if not self.parent:
            return False
        return self.parent.has_ancestor_of_type(sequence_type, include_self=True)

    def lift_over_to_first_ancestor_of_type(self, sequence_type: Union[str, SequenceType]) -> "Location":
        """Returns a new Location representing the liftover of this Location to its closest ancestor sequence (parent,
        parent of parent, etc.) which has the given sequence type. If the immediate parent has the given type,
        returns this Location. Raises NoSuchAncestorException if no ancestor with the given type exists."""
        try:
            self.first_ancestor_of_type(sequence_type)
        except NoSuchAncestorException:
            raise NoSuchAncestorException("Location has no ancestor of type {}".format(sequence_type))
        if self.parent_type == sequence_type:
            return self
        lifted_to_grandparent = self.parent.lift_child_location_to_parent()
        return lifted_to_grandparent.lift_over_to_first_ancestor_of_type(sequence_type)

    def has_ancestor_sequence(self, sequence: Sequence) -> bool:
        """Returns True iff this Location has some ancestor (parent, parent of parent, etc.) whose sequence
        attribute is equal to the given sequence"""
        if not self.parent:
            return False
        return self.parent.has_ancestor_sequence(sequence, include_self=True)

    def lift_over_to_sequence(self, sequence: Sequence) -> "Location":
        """Returns a new Location representing the liftover of this Location to the given sequence. The given
        sequence must be equal to the sequence attribute of some Parent in the ancestor hierarchy of this
        Location; otherwise, raises NoSuchAncestorException."""
        if not self.is_contiguous:
            raise ValueError("Location must be contiguous")
        if not self.has_ancestor_sequence(sequence):
            raise NoSuchAncestorException(
                "\nLocation:\n{}\nDoes not have ancestor:\n{}".format(str(self), sequence.summary())
            )
        if self.parent.sequence and self.parent.sequence == sequence:
            return self
        lifted_to_grandparent = self.parent.lift_child_location_to_parent()
        return lifted_to_grandparent.lift_over_to_sequence(sequence)

    def contains(
        self,
        other: "Location",
        match_strand: bool = False,
        full_span: bool = False,
        strict_parent_compare: bool = False,
    ) -> bool:
        """Returns True iff this location contains the other. If ``full_span`` is ``True``, the full span of
        both locations are compared.

        Parameters
        ----------
        other
            Other location
        match_strand
            If set to True, automatically return EmptyLocation() if other Location has a different Strand than
            this Location
        full_span
            If set to True, compare the full span of this Location to the full span of the other Location.
        strict_parent_compare
            Raise MismatchedParentException if parents do not match

        """
        if strict_parent_compare:
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        if not self.has_overlap(other, match_strand, full_span):
            return False
        other_to_compare = other.reset_parent(None)
        if full_span is False:
            return len(
                self.reset_parent(None).intersection(other_to_compare, match_strand=match_strand, full_span=full_span)
            ) == len(other_to_compare)
        else:
            return self.reset_parent(None)._full_span_interval.contains(other_to_compare._full_span_interval)
