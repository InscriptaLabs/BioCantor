__version__ = "0.5.0"

from abc import ABC, abstractmethod
from enum import Enum
from typing import TypeVar, List, Iterator, Union, Optional

from Bio.SeqFeature import FeatureLocation, CompoundLocation

Strand = TypeVar("Strand")
Alphabet = TypeVar("Alphabet")


class SequenceType(str, Enum):
    CHROMOSOME = "chromosome"
    SEQUENCE_CHUNK = "sequence_chunk"


class DistanceType(Enum):
    INNER = "inner"
    OUTER = "outer"
    STARTS = "starts"
    ENDS = "ends"


class AbstractLocation(ABC):
    """Shared AbstractLocation base class simplifies imports for type checking"""

    # The 0-based start position of this Location on its parent
    start: int

    # The 0-based exclusive end position of this Location on its parent
    end: int

    # The strand of this Location with respect to its parent
    strand: Strand

    # The parent of this Location
    parent: "AbstractParent"

    # The length (number of positions) of this Location. For subclasses representing discontiguous locations,
    # regions between blocks are not considered
    length: int

    def __len__(self):
        """Returns the length (number of positions) of this Location. For subclasses representing discontiguous
        locations, regions between blocks are not considered."""
        return self.length

    @abstractmethod
    def __str__(self):
        """Returns a human readable string representation of this Location"""

    @abstractmethod
    def __eq__(self, other):
        """Returns True iff this Location is equal to other object"""

    @abstractmethod
    def __hash__(self):
        """Returns a hash code satisfying location1 == location2 => hash(location1) == hash(location2)"""

    @abstractmethod
    def __repr__(self):
        """Returns the 'official' string representation of this Location"""

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
    def is_contiguous(self) -> bool:
        """Returns True iff this Location is fully contiguous within its parent"""

    @property
    @abstractmethod
    def is_empty(self) -> bool:
        """Returns True iff this Location is empty"""

    @property
    @abstractmethod
    def blocks(self) -> List["AbstractLocation"]:
        """Returns list of contiguous blocks comprising this Location"""

    @property
    @abstractmethod
    def is_overlapping(self) -> bool:
        """Returns True if this interval contains overlaps; always False for SingleInterval"""

    @property
    @abstractmethod
    def _full_span_interval(self) -> "AbstractLocation":
        """Returns the full span of this interval; is trivial for a SingleInterval and EmptyLocation"""

    @abstractmethod
    def scan_blocks(self) -> Iterator["AbstractLocation"]:
        """Returns an iterator over blocks in order relative to strand of this Location"""

    @property
    @abstractmethod
    def num_blocks(self) -> int:
        """Returns number of contiguous blocks comprising this Location"""

    @abstractmethod
    def optimize_blocks(self) -> "AbstractLocation":
        """Returns a new Location covering the same positions but with blocks optimized.
        For example, empty blocks may be removed or adjacent blocks may be combined if applicable."""

    @abstractmethod
    def gap_list(self) -> List["AbstractLocation"]:
        """Returns list of contiguous regions comprising the space between blocks of this Location. List is
        ordered relative to strand of this Location."""

    @abstractmethod
    def gaps_location(self) -> "AbstractLocation":
        """Returns a Location representing the space between blocks of this Location."""

    @abstractmethod
    def extract_sequence(self) -> "AbstractSequence":
        """Extracts the sequence of this Location from the parent.
        Concrete implementations should raise ValueError if no parent exists."""

    @abstractmethod
    def parent_to_relative_pos(self, parent_pos: int) -> int:
        """Converts a position on the parent to a position relative to this Location.
        Concrete implementations should raise ValueError if the given position does not overlap this Location."""

    @abstractmethod
    def relative_to_parent_pos(self, relative_pos: int) -> int:
        """Converts a position relative to this Location to a position on the parent"""

    @abstractmethod
    def parent_to_relative_location(
        self, parent_location: "AbstractLocation", optimize_blocks: bool = True
    ) -> "AbstractLocation":
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

    def location_relative_to(self, other: "AbstractLocation", optimize_blocks: bool = True) -> "AbstractLocation":
        """Converts this Location to a Location relative to another Location. The Locations must overlap.
        The returned value represents the relative location of the overlap within the other Location.

        If ``optimize_blocks`` is ``True``, the resulting Location will not have any adjacent or overlapping
        intervals. This is often desirable, because the output of this function can have weird coordinates
        when the locations are overlapping or adjacent. However, there are some cases where it is desirable
        to retain the original block structure. One such example are CDS where adjacent blocks or overlapping
        blocks are used to model frameshifts or indels.
        """

    @abstractmethod
    def _location_relative_to(self, other: "AbstractLocation", optimize_blocks: bool = True) -> "AbstractLocation":
        raise NotImplementedError

    @abstractmethod
    def relative_interval_to_parent_location(
        self, relative_start: int, relative_end: int, relative_strand: Strand
    ) -> "AbstractLocation":
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

    @abstractmethod
    def scan_windows(self, window_size: int, step_size: int, start_pos: int = 0) -> Iterator["AbstractLocation"]:
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

    @abstractmethod
    def has_overlap(
        self,
        other: "AbstractLocation",
        match_strand: bool = False,
        full_span: bool = False,
        strict_parent_compare: bool = False,
    ) -> bool:
        """Returns True iff this Location shares at least one position with the given Location.
        For subclasses representing discontiguous locations, regions between blocks are not considered.

        Parameters
        ----------
        other
            Other Location
        match_strand
            If set to True, automatically return False if given interval Strand does not match this Location's Strand
        full_span
            If set to True, compare the full span of this Location to the full span of the other Location.
        strict_parent_compare
            Raise MismatchedParentException if parents do not match

        Returns
        -------
        True if there is any overlap, False otherwise
        """

    @abstractmethod
    def reverse(self) -> "AbstractLocation":
        """Returns a new Location corresponding to this Location with the same start and stop, with
        strand and structure reversed"""

    @abstractmethod
    def reverse_strand(self) -> "AbstractLocation":
        """Returns a new Location corresponding to this Location with the strand reversed"""

    @abstractmethod
    def reset_strand(self, new_strand: Strand) -> "AbstractLocation":
        """Returns a new Location corresponding to this Location with the given strand"""

    @abstractmethod
    def reset_parent(self, new_parent: Optional["AbstractParent"]) -> "AbstractLocation":
        """Returns a new Location corresponding to this Location with positions unchanged and pointing
        to a new parent"""

    @abstractmethod
    def shift_position(self, shift: int) -> "AbstractLocation":
        """Returns a new Location corresponding to this location shifted by the given distance"""

    @abstractmethod
    def distance_to(self, other: "AbstractLocation", distance_type: DistanceType = DistanceType.INNER) -> int:
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
    def merge_overlapping(self) -> "AbstractLocation":
        """Merges overlapping windows"""

    @abstractmethod
    def to_biopython(self) -> Union[FeatureLocation, CompoundLocation]:
        """Returns a BioPython interval type; since they do not have a shared base class, we need a union"""

    @abstractmethod
    def first_ancestor_of_type(self, sequence_type: Union[str, SequenceType]) -> "AbstractParent":
        """Returns the Parent object representing the closest ancestor (parent, parent of parent, etc.)
        of this location which has the given sequence type. Raises NoSuchAncestorException if no ancestor with
        the given type exists."""

    @abstractmethod
    def has_ancestor_of_type(self, sequence_type: Union[str, SequenceType]) -> bool:
        """Returns True if some ancestor (parent, parent of parent, etc.) of of this location has the given sequence
        type, or False otherwise."""

    @abstractmethod
    def lift_over_to_first_ancestor_of_type(self, sequence_type: Union[str, SequenceType]) -> "AbstractLocation":
        """Returns a new Location representing the liftover of this Location to its closest ancestor sequence (parent,
        parent of parent, etc.) which has the given sequence type. If the immediate parent has the given type,
        returns this Location. Raises NoSuchAncestorException if no ancestor with the given type exists."""

    def has_ancestor_sequence(self, sequence: "AbstractSequence") -> bool:
        """Returns True iff this Location has some ancestor (parent, parent of parent, etc.) whose sequence
        attribute is equal to the given sequence"""

    def lift_over_to_sequence(self, sequence: "AbstractSequence") -> "AbstractLocation":
        """Returns a new Location representing the liftover of this Location to the given sequence. The given
        sequence must be equal to the sequence attribute of some Parent in the ancestor hierarchy of this
        Location; otherwise, raises NoSuchAncestorException."""

    @abstractmethod
    def intersection(
        self,
        other: "AbstractLocation",
        match_strand: bool = True,
        full_span: bool = False,
        strict_parent_compare: bool = False,
    ) -> "AbstractLocation":
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
        full_span
            If set to True, compare the full span of this Location to the full span of the other Location.
        strict_parent_compare
            Raise MismatchedParentException if parents do not match

        """

    @abstractmethod
    def union(self, other: "AbstractLocation") -> "AbstractLocation":
        """Returns a new Location representing the union of this Location with the other Location. This operation
        is commutative. Raises exception if locations cannot be combined."""

    @abstractmethod
    def union_preserve_overlaps(self, other: "AbstractLocation") -> "AbstractLocation":
        """Returns a new Location representing the union of this Location with the other Location, retaining
        overlapping blocks where applicable. This operation is commutative. Raises exception if locations cannot
        be combined."""

    @abstractmethod
    def minus(
        self, other: "AbstractLocation", match_strand: bool = True, strict_parent_compare: bool = False
    ) -> "AbstractLocation":
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
        strict_parent_compare
            Raise MismatchedParentException if parents do not match

        """

    @abstractmethod
    def extend_absolute(self, extend_start: int, extend_end: int) -> "AbstractLocation":
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
    def extend_relative(self, extend_upstream: int, extend_downstream: int) -> "AbstractLocation":
        """Returns a new Location extended upstream and downstream relative to this Location's Strand.

        Parameters
        ----------
        extend_upstream
            Non-negative integer: amount to extend upstream relative to Strand
        extend_downstream
            Non-negative integer: amount to extend downstream relative to Strand
        """

    @abstractmethod
    def contains(
        self,
        other: "AbstractLocation",
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


class AbstractSequence(ABC):
    """Shared AbstractSequence base class simplifies imports for type checking"""

    sequence_type: SequenceType
    _len: int
    sequence: str
    alphabet: Alphabet
    parent: Optional["AbstractParent"]
    id: Optional[str]

    def __len__(self):
        return self._len


class AbstractParent(ABC):
    """Shared AbstractParent base class simplifies imports for type checking"""

    id: Optional[str]
    sequence_type: Optional[SequenceType]
    sequence: Optional[AbstractSequence]
    location: Optional[AbstractLocation]
    parent: Optional["AbstractParent"]
    _strand: Optional[Strand]

    @abstractmethod
    def equals_except_location(self, other, require_same_sequence: bool = True):
        """Checks that this Parent is equal to another Parent, ignoring the associated Location members.

        By default also checks that any associated Sequence objects also match, but this can be toggled off.
        """

    @property
    @abstractmethod
    def strand(self) -> Optional[Strand]:
        """Returns the Strand of this Parent. If this Parent has no explicit Strand, but has a Location,
        that Location's Strand is returned."""

    @abstractmethod
    def strip_location_info(self) -> "AbstractParent":
        """Returns a new Parent object representing this Parent with information about child
        location removed"""

    @abstractmethod
    def first_ancestor_of_type(
        self, sequence_type: Union[str, SequenceType], include_self: bool = True
    ) -> "AbstractParent":
        """Returns the Parent object representing the closest ancestor (parent, parent of parent, etc.)
        of this Parent which has the given sequence type. If include_self is True and this Parent
        has the given type, returns this object. Raises NoSuchAncestorException if no ancestor with the given
        type exists.

        Parameters
        ----------
        sequence_type: str
            Sequence type
        include_self:
            Include this sequence as a candidate
        """

    @abstractmethod
    def has_ancestor_of_type(self, sequence_type: Union[str, SequenceType], include_self: bool = True) -> bool:
        """Returns True if some ancestor (parent, parent of parent, etc.) of this Parent has the given sequence type,
        or False otherwise. If include_self is True and this Parent has the given type, returns True.

        Parameters
        ----------
        sequence_type: str
            Sequence type
        include_self:
            Include this sequence as a candidate
        """

    @abstractmethod
    def lift_child_location_to_parent(self):
        """Lifts location of child object on this parent to the parent of this parent.
        Raises ValueError if any required data is missing (child location or location of this parent
        on its parent).

        Returns
        -------
        Location
            Child object location lifted to the parent of this parent
        """

    @abstractmethod
    def reset_location(self, location) -> "AbstractParent":
        """Returns a new Parent object with child location set to the given location"""

    @abstractmethod
    def has_ancestor_sequence(self, sequence, include_self: bool = True) -> bool:
        """Returns True iff this Parent has some ancestor (parent, parent of parent, etc.) whose sequence
        attribute is equal to the given sequence. If include_self is True and this Parent has sequence
        attribute equal to the given sequence, returns True."""
