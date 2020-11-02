from functools import reduce
from typing import TypeVar, Optional

import inscripta.biocantor
from inscripta.biocantor.exc import NoSuchAncestorException
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.util.object_validation import ObjectValidation

Parent = TypeVar("Parent")
ParentInputType = TypeVar("ParentInputType")


class Parent:
    """
    Holds information about a parent of some object. Typically the child object should hold
    a reference to this parent.
    """

    def __init__(
        self,
        *,
        id: Optional[str] = None,
        sequence_type: Optional[str] = None,
        strand: Optional[Strand] = None,
        location: Optional = None,
        sequence: Optional = None,
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

        location_parent_id = location.parent_id if location else None
        sequence_id = sequence.id if sequence else None
        non_null_ids = set([x for x in [id, location_parent_id, sequence_id] if x is not None])
        if len(non_null_ids) > 1:
            raise ValueError(
                "ID, location parent ID, and sequence ID do not match: {}, {}, {}".format(
                    id, location_parent_id, sequence_id
                )
            )

        location_parent_type = location.parent_type if location else None
        sequence_seqtype = sequence.sequence_type if sequence else None
        non_null_types = set([x for x in [sequence_type, location_parent_type, sequence_seqtype] if x is not None])
        if len(non_null_types) > 1:
            raise ValueError(
                "Sequence type, location parent type, and sequence do not match: {}, {}, {}".format(
                    sequence_type, location_parent_type, sequence_seqtype
                )
            )

        if location:
            if strand and location.strand and strand != location.strand:
                raise ValueError("Strand does not match location: {} != {}".format(strand, location.strand))
            if sequence and location.end > len(sequence):
                raise ValueError(
                    "Location end ({}) is greater than sequence length ({})".format(location.end, len(sequence))
                )

        parent_obj = inscripta.biocantor.parent.make_parent(parent) if parent else None
        if sequence and parent_obj and parent_obj.sequence and len(sequence) > len(parent_obj.sequence):
            raise ValueError(
                "Parent ({}) is longer than parent of parent ({})".format(len(sequence), len(parent_obj.sequence))
            )

        if sequence and sequence.parent:
            if parent_obj:
                ObjectValidation.require_parents_equal_except_location(parent_obj, sequence.parent)
                self.parent = parent_obj
            else:
                self.parent = sequence.parent
        else:
            self.parent = parent_obj

        self._id = non_null_ids.pop() if non_null_ids else None
        self._sequence_type = non_null_types.pop() if non_null_types else None
        self._strand = strand
        self.location = location
        self.sequence = sequence

    def __eq__(self, other):
        if not self.equals_except_location(other):
            return False
        return self.location == other.location and self.strand == other.strand

    def equals_except_location(self, other):
        if type(other) is not Parent:
            return False
        if self.id != other.id:
            return False
        if self.sequence_type != other.sequence_type:
            return False
        if self.sequence != other.sequence:
            return False
        if self.parent and other.parent and self.parent != other.parent:
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
    def id(self):
        if self._id:
            return self._id
        if self.sequence:
            return self.sequence.id
        return None

    @property
    def sequence_type(self):
        if self._sequence_type:
            return self._sequence_type
        if self.sequence:
            return self.sequence.sequence_type
        return None

    @property
    def strand(self):
        if self._strand:
            return self._strand
        if self.location:
            return self.location.strand
        return None

    def strip_location_info(self) -> Parent:
        """Returns a new Parent object representing this Parent with information about child
        location removed"""
        return Parent(
            id=self.id,
            sequence_type=self.sequence_type,
            sequence=self.sequence,
            parent=self.parent,
        )

    def first_ancestor_of_type(self, sequence_type, include_self: bool = True) -> Parent:
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

    def has_ancestor_of_type(self, sequence_type, include_self: bool = True) -> bool:
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
        lifted_blocks = [
            self.parent.location.relative_interval_to_parent_location(block.start, block.end, block.strand)
            for block in self.location.blocks
        ]
        lifted_blocks_union = reduce(lambda location1, location2: location1.union(location2), lifted_blocks)
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
