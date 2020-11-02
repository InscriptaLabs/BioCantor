from enum import Enum
from typing import TypeVar, Optional, Union

from Bio.Seq import Seq
from inscripta.biocantor.exc import AlphabetError, NoSuchAncestorException, EmptySequenceFastaError
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import Parent, make_parent
from inscripta.biocantor.sequence.alphabet import (
    Alphabet,
    ALPHABET_TO_NUCLEOTIDE_COMPLEMENT,
)

Location = TypeVar("Location")
Sequence = TypeVar("Sequence")
ParentInputType = TypeVar("ParentInputType")


class Sequence:
    """A sequence with an alphabet"""

    sequence: Seq

    def __init__(
        self,
        data: str,
        alphabet: Alphabet,
        id: Optional[str] = None,
        type: Optional[str] = None,
        parent: Optional[ParentInputType] = None,
        validate_alphabet: bool = True,
    ):
        """
            Parameters
        ----------
            data
                The contents of the sequence
            alphabet
                Alphabet
            id
                Sequence name
            type
                Type of sequence
            parent
                Parent
            validate_alphabet
                Whether to validate this sequence against its alphabet
        """
        self.sequence = Seq(data)
        self.alphabet = alphabet
        self.id = id
        self.sequence_type = type
        self.parent = make_parent(parent) if parent else None
        self._len = len(self.sequence)
        if self.parent and self.parent.location and len(self.parent.location) != len(self):
            raise ValueError(
                "Sequence length ({}) does not equal parent location length ({})".format(
                    len(self), len(self.parent.location)
                )
            )
        if validate_alphabet:
            self._validate_alphabet()

    def __eq__(self, other):
        if type(other) is not Sequence:
            return False
        if self.id != other.id:
            return False
        if self.sequence_type != other.sequence_type:
            return False
        if self.alphabet != other.alphabet:
            return False
        if self.parent != other.parent:
            return False
        return self.sequence == other.sequence

    def __hash__(self):
        return hash((self.id, self.sequence_type, self.alphabet, self.parent, self.sequence))

    def __str__(self):
        """Returns the sequence data as a string"""
        return str(self.sequence)

    def __len__(self):
        return self._len

    def __getitem__(self, key: Union[int, slice]) -> Sequence:
        """Returns a slice of the current Sequence as a new Sequence object"""
        subseq = str(self.sequence[key])

        if self.parent is not None and self.parent.location is not None:
            if isinstance(key, slice):
                rel_start = key.start
                rel_end = key.stop
            else:
                rel_start = key
                rel_end = key + 1
            new_parent_location = self.parent.location.relative_interval_to_parent_location(
                relative_start=rel_start, relative_end=rel_end, relative_strand=Strand.PLUS
            )
            new_parent = self.parent.reset_location(new_parent_location)
        else:
            new_parent = self.parent

        return Sequence(subseq, self.alphabet, type=self.sequence_type, parent=new_parent)

    def __repr__(self):
        return "<{}>".format(self.summary())

    def summary(self) -> str:
        """Returns a short string summary of this Sequence"""
        if self.id:
            id = self.id
        else:
            if len(self) <= 20:
                id = "Sequence={}".format(str(self))
            else:
                id = "Sequence"
        return "{};\n  Alphabet={};\n  Length={};\n  Parent={};\n  Type={}".format(
            id, self.alphabet.name, len(self), repr(self.parent), self.sequence_type
        )

    def _validate_alphabet(self):
        """Raises AlphabetError if this Sequence does not conform to its alphabet"""
        Sequence.validate_alphabet(str(self), self.alphabet)

    @staticmethod
    def validate_alphabet(sequence, alphabet):
        if sequence.upper().strip(alphabet.value) != "":
            raise AlphabetError("Invalid sequence for alphabet {}".format(alphabet.name))

    @property
    def parent_id(self) -> Optional[str]:
        return self.parent.id if self.parent else None

    @property
    def is_empty(self) -> bool:
        """Is this a len 0 sequence?"""
        return self._len == 0

    @property
    def location_on_parent(self) -> Optional[Location]:
        """Location of this sequence relative to the parent"""
        return self.parent.location if self.parent else None

    @property
    def parent_strand(self) -> Optional[Strand]:
        """Strand of this sequence"""
        return self.parent.strand if self.parent else None

    @property
    def parent_type(self) -> Optional[str]:
        return self.parent.sequence_type if self.parent else None

    def reverse_complement(self, new_id: str = None, new_type: str = None) -> "Sequence":
        """Returns a new Sequence corresponding to the reverse complement of this Sequence.
        Location on parent, if it exists, is converted appropriately.

        Parameters
        ----------
        new_id
            ID for the returned Sequence. If no value is provided, None is used.
        new_type
            Sequence type for the returned Sequence. If no value is provided, None is used.
        """
        if not self.alphabet.is_nucleotide_alphabet():
            raise AlphabetError("Cannot reverse complement sequence with alphabet {}".format(self.alphabet))
        location = self.location_on_parent.reverse_strand() if self.location_on_parent else None
        strand = self.parent_strand.reverse() if self.parent_strand else None
        rc_map = ALPHABET_TO_NUCLEOTIDE_COMPLEMENT[self.alphabet]
        try:
            seq_data = "".join([rc_map[c] for c in str(self)[::-1]])
        except KeyError as e:
            raise AlphabetError("Character {} not found for alphabet {}".format(str(e), self.alphabet))
        rc_parent = Parent(strand=strand, location=location) if strand or location else None
        return Sequence(
            seq_data,
            self.alphabet,
            id=new_id,
            type=new_type,
            parent=rc_parent,
            validate_alphabet=False,
        )

    def append(self, other: Sequence, new_id: Optional[str] = None, data_only: bool = False) -> Sequence:
        """Returns a new Sequence consisting of other sequence appended to the end of this Sequence

        Parameters
        ----------
        other
            Other sequence
        new_id
            ID for the returned sequence
        data_only
            Only combine sequence data, ignoring and removing metadata such as sequence type and parent
            info. If False, all metadata for both sequences must be consistent, and the returned sequence
            will have appropriate combined metadata.
        """
        if self.alphabet != other.alphabet:
            raise ValueError("Sequences must have same alphabet: {} != {}".format(self.alphabet, other.alphabet))
        new_seq_data = "{}{}".format(str(self), str(other))
        if data_only:
            return Sequence(new_seq_data, self.alphabet, id=new_id)
        if self.sequence_type != other.sequence_type:
            raise ValueError("Sequences must have same type: {} != {}".format(self.sequence_type, other.sequence_type))
        if self.parent:
            if not self.parent.equals_except_location(other.parent):
                raise ValueError(
                    "Sequences must have same parent (except location on parent):\n{}\n  !=\n{}".format(
                        repr(self.parent), repr(other.parent)
                    )
                )
            if self.parent.strand == Strand.UNSTRANDED or self.parent.strand != other.parent.strand:
                raise ValueError("Invalid strands on parent: {} != {}".format(self.parent.strand, other.parent.strand))
            if self.parent.location and other.parent.location:
                if self.parent.strand == Strand.PLUS and self.parent.location.end > other.parent.location.start:
                    raise ValueError("Sequence on plus strand of parent must be to the left of appended sequence")
                if self.parent.strand == Strand.MINUS and self.parent.location.start < other.parent.location.end:
                    raise ValueError("Sequence on minus strand of parent must be to the right of appended sequence")
                new_location = self.parent.location.union(other.parent.location)
            else:
                new_location = None
            new_parent = self.parent.strip_location_info().reset_location(new_location)
        else:
            new_parent = None
        return Sequence(
            new_seq_data,
            self.alphabet,
            id=new_id,
            type=self.sequence_type,
            parent=new_parent,
            validate_alphabet=False,
        )

    def first_ancestor_of_type(self, sequence_type: str, include_self: bool = True) -> Parent:
        """Returns the Parent object representing the closest ancestor (parent, parent of parent, etc.)
        of this sequence which has the given sequence type. If include_self is True and this sequence has the
        given type, returns a new Parent object representing this sequence. Raises NoSuchAncestorException
        if no ancestor with the given type exists.

        Parameters
        ----------
        sequence_type:
            Sequence type
        include_self:
            Include this sequence as a candidate
        """
        if include_self and self.sequence_type == sequence_type:
            return Parent(sequence=self)
        if self.parent:
            return self.parent.first_ancestor_of_type(sequence_type, True)
        raise NoSuchAncestorException

    def has_ancestor_of_type(self, sequence_type: str, include_self: bool = True) -> bool:
        """Returns True if some ancestor (parent, parent of parent, etc.) of this sequence which has the given sequence
        type, or False otherwise. If include_self is True and this sequence has the given type, returns True.

        Parameters
        ----------
        sequence_type:
            Sequence type
        include_self:
            Include this sequence as a candidate
        """
        if include_self and self.sequence_type == sequence_type:
            return True
        if self.parent:
            return self.parent.has_ancestor_of_type(sequence_type, include_self=True)
        return False

    def to_fasta(self, num_chars: Optional[int] = 60) -> str:
        """Returns a FASTA-formatted string for this sequence. These are line-broken every num_chars.

        Parameters
        ----------
        num_chars:
            Number of characters per line. Defaults to 60, which is the same as BioPython.
        """
        if self.is_empty:
            raise EmptySequenceFastaError("Cannot write FASTA for empty Sequence")

        r = [f">{self.id}"]
        for i in range(0, self._len, num_chars):
            r.append(str(self)[i : i + num_chars])
        return "\n".join(r)
