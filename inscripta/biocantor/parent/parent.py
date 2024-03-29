from functools import reduce, lru_cache
from typing import TypeVar, Optional, Union, Iterable

# this base import is required in order to avoid a circular import; this module needs make_parent()
# it is not possible to put make_parent() within this module, because it references Parent, and so it becomes an
# internally circular reference
import inscripta.biocantor
from inscripta.biocantor import AbstractParent, AbstractLocation, AbstractSequence, SequenceType
from inscripta.biocantor.exc import (
    NoSuchAncestorException,
    LocationException,
    InvalidStrandException,
    ParentException,
    InvalidPositionException,
)
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.util.object_validation import ObjectValidation

Parent = TypeVar("Parent")
ParentInputType = TypeVar("ParentInputType")

# 1000 seems reasonable for Parent caches
# When parsing annotation files, the number of Parent objects built will likely be the # of chromosomes in the genome
# When sequence chunks are used, the number of Parents will equal the number of distinct chunks built
# In testing, a cache size of 1000 was more performant than 5000
PARENT_CACHE_SIZE = 1000


@lru_cache(maxsize=PARENT_CACHE_SIZE)
def _unique_value_or_none(values: Iterable[Optional[str]]) -> Optional[str]:
    """Checks if a set of values contains more than one distinct non-null value. If so, raises ValueError.
    Otherwise, returns the single unique non-null value (if there is one) or None if all values are None."""
    values = {x for x in values if x is not None}
    if len(values) == 1:
        return values.pop()
    elif len(values) == 0:
        return None
    else:
        raise ParentException(f"Multiple distinct non-null values were provided: {values}")


@lru_cache(maxsize=PARENT_CACHE_SIZE)
class Parent(AbstractParent):
    """
    Holds information about a parent of some object. Typically the child object should hold
    a reference to this parent.
    """

    __slots__ = []

    def __init__(
        self,
        *,
        id: Optional[str] = None,
        sequence_type: Optional[Union[SequenceType, str]] = None,
        strand: Optional[Strand] = None,
        location: Optional[AbstractLocation] = None,
        sequence: Optional[AbstractSequence] = None,
        parent: Optional[ParentInputType] = None,
    ):
        """
        Parameters
        ----------
        id
            Parent ID
        sequence_type: string
            Parent type
        strand: Strand
            Strand of child object with respect to this parent
        location: Location
            Location of child object with respect to this parent
        sequence: Sequence
            Sequence of this parent
        parent: Parent
            Parent of this parent
        """

        location_parent_id = location.parent_id if location is not None else None
        sequence_id = sequence.id if sequence is not None else None
        parent_id = _unique_value_or_none((id, location_parent_id, sequence_id))

        location_parent_type = location.parent_type if location is not None else None
        sequence_seqtype = sequence.sequence_type if sequence is not None else None
        seq_type = _unique_value_or_none((sequence_type, location_parent_type, sequence_seqtype))

        if location is not None:
            if strand and location.strand and strand is not location.strand:
                raise InvalidStrandException("Strand does not match location: {} != {}".format(strand, location.strand))
            if sequence is not None and location.end > len(sequence):
                raise InvalidPositionException(
                    "Location end ({}) is greater than sequence length ({})".format(location.end, len(sequence))
                )

        parent_obj = inscripta.biocantor.parent.make_parent(parent) if parent else None
        if (
            sequence is not None
            and parent_obj
            and parent_obj.sequence is not None
            and len(sequence) > len(parent_obj.sequence)
        ):
            raise LocationException(
                "Parent ({}) is longer than parent of parent ({})".format(len(sequence), len(parent_obj.sequence))
            )

        if sequence is not None and sequence.parent is not None:
            if parent_obj:
                ObjectValidation.require_parents_equal_except_location(parent_obj, sequence.parent)
                self.parent = parent_obj
            else:
                self.parent = sequence.parent
        else:
            self.parent = parent_obj

        self.id = parent_id
        self.sequence_type = seq_type
        self._strand = strand
        self.location = location
        self.sequence = sequence

        self._strand_property = None

    def __eq__(self, other):
        if not self.equals_except_location(other):
            return False
        return self.location == other.location and self.strand is other.strand

    def equals_except_location(self, other, require_same_sequence: bool = True):
        """Checks that this Parent is equal to another Parent, ignoring the associated Location members.

        By default also checks that any associated Sequence objects also match, but this can be toggled off.
        """
        # this checks the object under the lru_cache hood
        if type(other) is not Parent.__wrapped__:
            return False
        if self.id != other.id:
            return False
        if self.sequence_type != other.sequence_type:
            return False
        if self.parent and other.parent and self.parent != other.parent:
            return False
        if require_same_sequence and self.sequence != other.sequence:
            return False
        return True

    def __hash__(self):
        return hash(
            (
                self.id,
                self.sequence_type,
                self.strand,
                self.location,
                self.sequence,
                self.parent,
            )
        )

    def __repr__(self):
        return "<Parent: id={}, type={}, strand={}, location={}, sequence={}, parent={}>".format(
            self.id,
            self.sequence_type,
            self.strand,
            repr(self.location),
            repr(self.sequence),
            repr(self.parent),
        )

    @property
    def strand(self) -> Optional[Strand]:
        """Returns the Strand of this Parent. If this Parent has no explicit Strand, but has a Location,
        that Location's Strand is returned."""
        if self._strand_property is None:
            if self._strand:
                self._strand_property = self._strand
            if self.location:
                self._strand_property = self.location.strand
        return self._strand_property

    def strip_location_info(self) -> Parent:
        """Returns a new Parent object representing this Parent with information about child
        location removed"""
        return Parent(
            id=self.id,
            sequence_type=self.sequence_type,
            sequence=self.sequence,
            parent=self.parent,
        )

    def first_ancestor_of_type(self, sequence_type: Union[str, SequenceType], include_self: bool = True) -> Parent:
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
        if include_self and self.sequence_type == sequence_type:
            return self
        if self.parent:
            return self.parent.first_ancestor_of_type(sequence_type)
        raise NoSuchAncestorException

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
        if include_self and self.sequence_type == sequence_type:
            return True
        if self.parent:
            return self.parent.has_ancestor_of_type(sequence_type, include_self=True)
        return False

    def lift_child_location_to_parent(self):
        """Lifts location of child object on this parent to the parent of this parent.
        Raises ValueError if any required data is missing (child location or location of this parent
        on its parent).

        Returns
        -------
        Location
            Child object location lifted to the parent of this parent
        """
        ObjectValidation.require_parent_has_location(self)
        ObjectValidation.require_parent_has_parent_with_location(self)
        lifted_blocks = (
            self.parent.location.relative_interval_to_parent_location(block.start, block.end, block.strand)
            for block in self.location.blocks
        )
        lifted_blocks_union = reduce(
            lambda location1, location2: location1.union_preserve_overlaps(location2), lifted_blocks
        )
        location_with_parent = lifted_blocks_union.reset_parent(self.parent.strip_location_info())
        return location_with_parent

    def reset_location(self, location) -> "Parent":
        """Returns a new Parent object with child location set to the given location"""
        strand = location.strand if location else None
        return Parent(
            id=self.id,
            sequence_type=self.sequence_type,
            strand=strand,
            location=location,
            sequence=self.sequence,
            parent=self.parent,
        )

    def has_ancestor_sequence(self, sequence, include_self: bool = True) -> bool:
        """Returns True iff this Parent has some ancestor (parent, parent of parent, etc.) whose sequence
        attribute is equal to the given sequence. If include_self is True and this Parent has sequence
        attribute equal to the given sequence, returns True."""
        if include_self and self.sequence and self.sequence == sequence:
            return True
        if not self.parent:
            return False
        return self.parent.has_ancestor_sequence(sequence, include_self=True)
