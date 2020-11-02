from abc import ABC, abstractmethod
from typing import Iterator, List, Union

from Bio.SeqFeature import FeatureLocation, CompoundLocation

from inscripta.biocantor.exc import NoSuchAncestorException
from inscripta.biocantor.location.distance import DistanceType
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import Parent
from inscripta.biocantor.sequence import Sequence
from inscripta.biocantor.util.object_validation import ObjectValidation


class Location(ABC):
    """Abstract location with respect to a coordinate system"""

    @abstractmethod
    def __len__(self):
        """Returns the length (number of positions) of this Location. For subclasses representing discontiguous
        locations, regions between blocks are not considered."""

    @abstractmethod
    def __str__(self):
        """Returns a human readable string representation of this Location"""

    @abstractmethod
    def __eq__(self, other):
        """Returns True iff this Location is equal to other object"""

    @abstractmethod
    def __repr__(self):
        """Returns the 'official' string representation of this Location"""

    @property
    @abstractmethod
    def parent(self) -> Parent:
        """Returns the parent of this Location"""

    @property
    def parent_id(self) -> str:
        """Returns the parent ID"""
        return self.parent.id if self.parent else None

    @property
    def parent_type(self) -> str:
        """Returns the sequence type of the parent"""
        return self.parent.sequence_type if self.parent else None

    @property
    @abstractmethod
    def strand(self) -> Strand:
        """Returns the strand of this Location with respect to its parent"""

    @property
    @abstractmethod
    def start(self) -> int:
        """Returns the 0-based start position of this Location on its parent"""

    @property
    @abstractmethod
    def end(self) -> int:
        """Returns the 0-based exclusive end position of this Location on its parent"""

    @property
    @abstractmethod
    def is_contiguous(self) -> bool:
        """Returns True iff this Location is fully contiguous within its parent"""

    @property
    @abstractmethod
    def is_empty(self) -> bool:
        """Returns True iff this Location is empty"""

    @property
    @abstractmethod
    def blocks(self) -> List["Location"]:
        """Returns list of contiguous blocks comprising this Location"""

    @property
    @abstractmethod
    def is_overlapping(self) -> bool:
        """Returns True if this interval contains overlaps; always False for SingleInterval"""

    @abstractmethod
    def scan_blocks(self) -> Iterator["Location"]:
        """Returns an iterator over blocks in order relative to strand of this Location"""

    @property
    @abstractmethod
    def num_blocks(self) -> int:
        """Returns number of contiguous blocks comprising this Location"""

    @abstractmethod
    def optimize_blocks(self) -> "Location":
        """Returns a new Location covering the same positions but with blocks optimized.
        For example, empty blocks may be removed or adjacent blocks may be combined if applicable."""

    @abstractmethod
    def gap_list(self) -> List["Location"]:
        """Returns list of contiguous regions comprising the space between blocks of this Location. List is
        ordered relative to strand of this Location."""

    @abstractmethod
    def gaps_location(self) -> "Location":
        """Returns a Location representing the space between blocks of this Location."""

    @abstractmethod
    def extract_sequence(self) -> Sequence:
        """Extracts the sequence of this Location from the parent.
        Concrete implementations should raise ValueError if no parent exists."""

    @abstractmethod
    def parent_to_relative_pos(self, parent_pos: int) -> int:
        """Converts a position on the parent to a position relative to this Location.
        Concrete implementations should raise ValueError if the given position does not overlap this Location."""

    @abstractmethod
    def relative_to_parent_pos(self, relative_pos: int) -> int:
        """Converts a position relative to this Location to a position on the parent"""

    def parent_to_relative_location(self, parent_location: "Location") -> "Location":
        """Converts a Location on the parent to a Location relative to this Location.

        Parameters
        ----------
        parent_location
            Location with the same parent as this Location. Both parents can be None.

        Returns
        -------
        New Location relative to this Location.
        """
        return parent_location.location_relative_to(self)

    def location_relative_to(self, other: "Location") -> "Location":
        """Converts this Location to a Location relative to another Location. The Locations must overlap.
        The returned value represents the relative location of the overlap within the other Location."""
        if other.is_empty:
            return other
        if self.parent or other.parent:
            if not self.parent and other.parent:
                raise ValueError(
                    "Parents must be both null or both non-null:\n{}\n  !=\n{}".format(self.parent, other.parent)
                )
            ObjectValidation.require_parents_equal_except_location(self.parent, other.parent)
        ObjectValidation.require_locations_overlap(self, other)
        return self._location_relative_to(other)

    @abstractmethod
    def _location_relative_to(self, other: "Location") -> "Location":
        raise NotImplementedError

    @abstractmethod
    def relative_interval_to_parent_location(
        self, relative_start: int, relative_end: int, relative_strand: Strand
    ) -> "Location":
        """Converts an interval relative to this Location to a Location on the parent

        Parameters
        ----------
        relative_start
            0-based start position of interval relative to this Location
        relative_end
            0-based exclusive end position of interval relative to this Location
        relative_strand
            Strand of interval relative to the strand of this Location. If the strand of interval is on the SAME
            strand as the strand of this location, relative_strand is PLUS. If the strand interval is on the OPPOSITE
            strand, relative_strand is MINUS.

        Returns
        -------
        New Location on the parent with the parent as parent
        """

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

    @abstractmethod
    def has_overlap(self, other: "Location", match_strand: bool = False) -> bool:
        """Returns True iff this Location shares at least one position with the given Location.
        For subclasses representing discontiguous locations, regions between blocks are not considered.

        Parameters
        ----------
        other
            Other Location
        match_strand
            If set to True, automatically return False if given interval Strand does not match this Location's Strand

        Returns
        -------
        True if there is any overlap, False otherwise
        """

    @abstractmethod
    def reverse(self) -> "Location":
        """Returns a new Location corresponding to this Location with the same start and stop, with
        strand and structure reversed"""

    @abstractmethod
    def reverse_strand(self) -> "Location":
        """Returns a new Location corresponding to this Location with the strand reversed"""

    @abstractmethod
    def reset_strand(self, new_strand: Strand) -> "Location":
        """Returns a new Location corresponding to this Location with the given strand"""

    @abstractmethod
    def reset_parent(self, new_parent: Parent) -> "Location":
        """Returns a new Location corresponding to this Location with positions unchanged and pointing
        to a new parent"""

    @abstractmethod
    def shift_position(self, shift: int) -> "Location":
        """Returns a new Location corresponding to this location shifted by the given distance"""

    @abstractmethod
    def distance_to(self, other: "Location", distance_type: DistanceType = DistanceType.INNER) -> int:
        """Returns the distance from this location to another location with the same parent.
        Return value is a non-negative integer and implementations must be commutative.

        Parameters
        ----------
        other
            Other location with same parent as this location
        distance_type
            Distance type
        """

    @abstractmethod
    def merge_overlapping(self) -> "Location":
        """Merges overlapping windows"""

    @abstractmethod
    def to_biopython(self) -> Union[FeatureLocation, CompoundLocation]:
        """Returns a BioPython interval type; since they do not have a shared base class, we need a union"""
        pass

    def first_ancestor_of_type(self, sequence_type: str) -> Parent:
        """Returns the Parent object representing the closest ancestor (parent, parent of parent, etc.)
        of this location which has the given sequence type. Raises NoSuchAncestorException if no ancestor with
        the given type exists."""
        if not self.parent:
            raise NoSuchAncestorException("Location has no parent")
        return self.parent.first_ancestor_of_type(sequence_type, include_self=True)

    def has_ancestor_of_type(self, sequence_type: str) -> bool:
        """Returns True if some ancestor (parent, parent of parent, etc.) of of this location has the given sequence
        type, or False otherwise."""
        if not self.parent:
            return False
        return self.parent.has_ancestor_of_type(sequence_type, include_self=True)

    def lift_over_to_first_ancestor_of_type(self, sequence_type: str) -> "Location":
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

    @abstractmethod
    def intersection(self, other: "Location", match_strand: bool = True) -> "Location":
        """Returns a new Location representing the intersection of this Location with the other Location.
        Returned Location, if nonempty, has the same Strand as this Location. This operation is commutative
        if match_strand is True.

        Parameters
        ----------
        other
            Other location
        match_strand
            If set to True, automatically return EmptyLocation() if other Location has a different Strand than
            this Location
        """

    @abstractmethod
    def union(self, other: "Location") -> "Location":
        """Returns a new Location representing the union of this Location with the other Location. This operation
        is commutative. Raises ValueError if locations cannot be combined."""

    @abstractmethod
    def minus(self, other: "Location", match_strand: bool = True) -> "Location":
        """Returns a new Location representing this Location minus its intersection with the other Location.
        Returned Location has the same Strand as this Location. If there is no intersection, returns this Location.
        This operation is not commutative.

        Parameters
        ----------
        other
            Other location
        match_strand
            If set to True, automatically return this Location if other Location has a different Strand than this
            Location
        """

    @abstractmethod
    def extend_absolute(self, extend_start: int, extend_end: int) -> "Location":
        """Returns a new Location representing this Location with start and end positions extended by the given
        values, ignoring Strand. Returned Location has same Strand as this Location.

        Parameters
        ----------
        extend_start
            Non-negative integer: amount to extend start
        extend_end
            Non-negative integer: amount to extend end
        """

    @abstractmethod
    def extend_relative(self, extend_upstream: int, extend_downstream: int) -> "Location":
        """Returns a new Location extended upstream and downstream relative to this Location's Strand.

        Parameters
        ----------
        extend_upstream
            Non-negative integer: amount to extend upstream relative to Strand
        extend_downstream
            Non-negative integer: amount to extend downstream relative to Strand
        """

    def contains(self, other: "Location", match_strand: bool = False):
        """Returns True iff this location contains the other"""
        if not self.has_overlap(other, match_strand):
            return False
        other_to_compare = other.reset_parent(None)
        return len(self.reset_parent(None).intersection(other_to_compare, match_strand=match_strand)) == len(
            other_to_compare
        )
