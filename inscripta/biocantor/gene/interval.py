"""
This module contains abstract base classes for interval types and interval collection types.
"""
from abc import ABC, abstractmethod
from enum import Enum
from typing import List, Union, Dict, Hashable, Set, Optional, Any, Iterable, Iterator, TypeVar, TYPE_CHECKING
from uuid import UUID

from methodtools import lru_cache

from inscripta.biocantor.exc import (
    ValidationException,
    NullSequenceException,
    NullParentException,
    NoSuchAncestorException,
    LocationOverlapException,
)
from inscripta.biocantor.io.bed import RGB, BED12
from inscripta.biocantor.io.gff3.rows import GFFRow
from inscripta.biocantor.location import Location, Strand
from inscripta.biocantor.location.location_impl import SingleInterval, CompoundInterval, EmptyLocation
from inscripta.biocantor.parent import Parent, SequenceType
from inscripta.biocantor.sequence import Sequence
from inscripta.biocantor.util.object_validation import ObjectValidation

# primitive data types possible as values of the list in a qualifiers dictionary
QualifierValue = TypeVar("QualifierValue", str, int, bool, float)

if TYPE_CHECKING:
    from inscripta.biocantor.gene.transcript import TranscriptInterval
    from inscripta.biocantor.gene.feature import FeatureInterval


class IntervalType(str, Enum):
    """This enum differentiates the three main types of Intervals -- Features, Transcripts and Variants"""

    FEATURE = "feature"
    TRANSCRIPT = "transcript"
    VARIANT = "variant"


class AbstractInterval(ABC):
    """This is a wrapper over :class:`~biocantor.location.Location` that adds metadata coordinate transformation
    QOL functions.

    All operations on coordinates are assumed to operate in chromosome–relative coordinates unless otherwise specified.
    All constructors use chromosome relative coordinates as well. If you want to operate on coordinate systems that
    are a subset of a chromosome, you must instantiate a Parent object that provides the coordinate relationship.

    A function to help you build these relationships can be found at :meth:`biocantor.io.parser.seq_chunk_to_parent()`.
    """

    _location: Location
    _identifiers: List[Union[str, UUID]]
    qualifiers: Dict[Hashable, Set[str]]  # all subclasses convert qualifier values to sets of strings
    guid: UUID
    sequence_guid: Optional[UUID] = None
    sequence_name: Optional[str] = None
    bin: int
    start: int
    end: int
    _parent_or_seq_chunk_parent: Optional[Parent] = None

    def __len__(self):
        return self.end - self.start

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        elif not self.to_dict() == other.to_dict():
            return False
        else:
            return self.chunk_relative_location == other.chunk_relative_location

    def __hash__(self):
        """Produces a hash, which is the GUID."""
        return hash((self.guid, self.chunk_relative_location))

    @abstractmethod
    def to_dict(self, chromosome_relative_coordinates: bool = True) -> Dict[str, Any]:
        """Dictionary to build Model representation. Defaults to always exporting in original chromosome
        relative coordinates, but this can be disabled to export in sequence-chunk relative coordinates.

        If you have exported to sequence-chunk relative coordinates, and then try to re-instantiate, the subsequent
        object will now consider these new coordinates to be the original chromosome coordinates, and the relationship
        back to the true coordinates will be lost.
        """

    def _parent_to_dict(self, chromosome_relative_coordinates: bool = True) -> Optional[Dict[str, Any]]:
        """
        Converts the ``_parent_or_seq_chunk_parent`` member of this Interval to a JSON-serializable representation.

        Raises:
            NotImplementedError: If chromosome_relative_coordinates is ``False`` and
                ``self._parent_or_seq_chunk_parent`` is not ``None``.
            NoSuchAncestorException: If ``self._parent_or_seq_chunk_parent`` is chunk-relative but lacks
                sequence information.
        """
        if not self._parent_or_seq_chunk_parent:
            return None
        elif chromosome_relative_coordinates is False:
            raise NotImplementedError("Cannot export parent to chunk relative coordinates")

        if self.chunk_relative_location.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            chunk_parent = self.chunk_relative_location.first_ancestor_of_type(SequenceType.SEQUENCE_CHUNK)
            sequence = chunk_parent.sequence
            if not sequence:
                raise NoSuchAncestorException("Chunk parents must have sequence")
            location = sequence.location_on_parent
            return {
                "seq": str(sequence),
                "sequence_name": self.chromosome_location.parent.id,
                "start": location.start,
                "end": location.end,
                "strand": location.strand.name,
                "alphabet": sequence.alphabet.name,
                "type": SequenceType.SEQUENCE_CHUNK.name,
            }
        elif self.chunk_relative_location.has_ancestor_of_type(SequenceType.CHROMOSOME):
            parent = self.chunk_relative_location.first_ancestor_of_type(SequenceType.CHROMOSOME)
            sequence = parent.sequence
            location = self.chromosome_location
            return {
                "seq": str(sequence) if sequence else None,
                "sequence_name": parent.id,
                "start": location.start,
                "end": location.end,
                "strand": location.strand.name,
                "alphabet": sequence.alphabet.name if sequence else None,
                "type": SequenceType.CHROMOSOME.name,
            }
        else:
            parent = self._parent_or_seq_chunk_parent
            seq_str = None
            start = None
            end = None
            strand = None
            alphabet = None
            sequence_type = None
            if parent.location:
                start = parent.location.start
                end = parent.location.end
                strand = parent.location.strand.name
            if parent.sequence:
                sequence_type = parent.sequence_type
                seq_str = str(parent.sequence)
                alphabet = parent.sequence.alphabet.name
            return {
                "seq": seq_str,
                "sequence_name": parent.id,
                "start": start,
                "end": end,
                "strand": strand,
                "alphabet": alphabet,
                "type": sequence_type,
            }

    @staticmethod
    @abstractmethod
    def from_dict(vals: Dict[str, Any], parent_or_seq_chunk_parent: Optional[Parent] = None) -> "AbstractInterval":
        """Build an interval from a dictionary representation"""

    @abstractmethod
    def to_gff(
        self,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterator[GFFRow]:
        """Writes a GFF format list of lists for this feature.

        Args:
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
        """

    @property
    def is_chunk_relative(self) -> bool:
        """Does this Interval object exist on a sequence chunk?"""
        return self.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK)

    @property
    def chunk_relative_size(self) -> int:
        return len(self.chunk_relative_location)

    @lru_cache(maxsize=1)
    @property
    def has_sequence(self) -> bool:
        """Returns true if this Interval has an associated sequence of any type"""
        try:
            ObjectValidation.require_location_has_parent_with_sequence(self.chunk_relative_location)
        except NullSequenceException:
            return False
        except NullParentException:
            return False
        else:
            return True

    @property
    @abstractmethod
    def id(self) -> str:
        """Returns the ID of this feature. Provides a shared API across genes/transcripts and features."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Returns the name of this feature. Provides a shared API across genes/transcripts and features."""

    @property
    def chunk_relative_start(self) -> int:
        """Returns chunk relative start position."""
        return self.chunk_relative_location.start

    @property
    def chunk_relative_end(self) -> int:
        """Returns chunk relative end position."""
        return self.chunk_relative_location.end

    @lru_cache(maxsize=1)
    @property
    def chromosome_location(self) -> Location:
        """Returns the Location of this in *chromosome coordinates*.

        If the coordinate system is unknown, this will return the same coordinate system as
        ``chunk_relative_location``, that is the true underlying ``_location`` member.

        This Location object will always have the full span of the Interval in chromosome coordinates,
        even if this feature exists in chunk relative coordinates. As a result of this, if this
        Interval was built on chunk relative coordinates, the sequence information will not be present.
        """
        if self._parent_or_seq_chunk_parent and self._parent_or_seq_chunk_parent.has_ancestor_of_type(
            SequenceType.CHROMOSOME
        ):
            parent = self._parent_or_seq_chunk_parent.first_ancestor_of_type(SequenceType.CHROMOSOME)
            return SingleInterval(self.start, self.end, Strand.PLUS, parent)
        else:
            return SingleInterval(self.start, self.end, Strand.PLUS)

    @lru_cache(maxsize=1)
    @property
    def _chunk_relative_bounded_chromosome_location(self) -> Location:
        """
        Returns the Location of this in *chromosome coordinates*.

        This function is different from ``chromosome_location`` in that it will return a Location bounded
        by the chunk relative location of this Interval, if it exists.

        This accessor is private because using it may lead to weird behavior. However, it is necessary
        for things like slicing CDSFrames in a chunk relative CDSInterval.
        """
        if self.chunk_relative_location.has_ancestor_of_type(SequenceType.CHROMOSOME):
            return self.lift_over_to_first_ancestor_of_type(SequenceType.CHROMOSOME)
        return self.chunk_relative_location

    @property
    def chunk_relative_location(self) -> Location:
        """Returns the Location of this in *chunk relative coordinates*"""
        return self._location

    @property
    def blocks(self) -> List[Location]:
        """Returns the blocks of this location"""
        return self.chromosome_location.blocks

    @property
    def num_blocks(self) -> int:
        """Returns the number of blocks of this location."""
        return self.chromosome_location.num_blocks

    @property
    def num_chunk_relative_blocks(self) -> int:
        """Returns the number of chunk-relative blocks of this location. Could be less than ``num_blocks``
        if this interval is a slice of the full length interval."""
        return self.chunk_relative_location.num_blocks

    @property
    def chunk_relative_blocks(self) -> List[Location]:
        """Returns the chunk relative blocks of this location"""
        return self.chunk_relative_location.blocks

    @property
    def strand(self) -> Strand:
        """Returns strand of location."""
        return self.chromosome_location.strand

    @property
    def chunk_relative_strand(self) -> Strand:
        """Returns strand of location."""
        return self.chunk_relative_location.strand

    @property
    def identifiers(self) -> Set[Union[str, UUID]]:
        """Returns the identifiers for this FeatureInterval, if they exist"""
        return {getattr(self, i) for i in self._identifiers if getattr(self, i) is not None}

    @property
    def identifiers_dict(self) -> Dict[str, Union[str, UUID]]:
        """Returns the identifiers and their keys for this FeatureInterval, if they exist"""
        return {key: getattr(self, key) for key in self._identifiers if getattr(self, key) is not None}

    @staticmethod
    def initialize_location(
        starts: List[int],
        ends: List[int],
        strand: Strand,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ) -> Location:
        """
        Initialize the :class:`Location` object for this interval.

        Args:
            starts: Start positions relative to the chromosome.
            ends: End positions relative to the chromosome.
            strand: Strand relative to the chromosome.
            parent_or_seq_chunk_parent: An optional parent, either as a full chromosome or as a sequence chunk.
        """
        if len(starts) != len(ends):
            raise ValidationException("Number of interval starts does not match number of interval ends")
        elif len(starts) == len(ends) == 1:
            location = SingleInterval(starts[0], ends[0], strand)
        else:
            location = CompoundInterval(starts, ends, strand)
        return AbstractInterval.liftover_location_to_seq_chunk_parent(location, parent_or_seq_chunk_parent)

    @staticmethod
    def liftover_location_to_seq_chunk_parent(
        location: Location,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ) -> Location:
        """
        BioCantor supports constructing any of the interval classes from a subset of the chromosome. In order to
        be able to set up the coordinate relationship and successfully pull down sequence, this function
        lifts the coordinates from the original annotation object on to this new coordinate system.

        .. code:: python

            parent_1_15 = Parent(
                sequence=Sequence(
                    genome2[1:15],
                    Alphabet.NT_EXTENDED_GAPPED,
                    type=SequenceType.SEQUENCE_CHUNK,
                    parent=Parent(
                        location=SingleInterval(1, 15, Strand.PLUS,
                                               parent=Parent(id="genome_1_15", sequence_type=SequenceType.CHROMOSOME))
                    ),
                )
            )

        Alternatively, if the sequence is coming straight from a file, it will be a :class:`Parent` with a
        :class:`Sequence` attached:

        .. code:: python

            parent = Parent(id="chr1", sequence=Sequence(genome, Alphabet.NT_STRICT, type=SequenceType.CHROMOSOME))

        This convenience function detects which kind of parent is given, and sets up the appropriate location.

        This function also handles the case where the ``location`` argument is already chunk-relative. If this is the
        case, the ``location`` object is first lifted back to its chromosomal coordinates, then lifted back down
        on to this new chunk.

        Args:
            location: A location object, likely produced by :meth:`initialize_location()`. Could also be the location
                of an existing AbstractInterval subclass, such as when the method
                ``liftover_interval_to_parent_or_seq_chunk_parent()`` is called.
            parent_or_seq_chunk_parent: An optional parent, either as a full chromosome or as a sequence chunk. If
                not provided, this function is a no-op.

        Returns:
            A :class:`Location` object.

        Raises:
            ValidationException: If ``parent_or_seq_chunk_parent`` has no ancestor of type ``chromosome`` or
                ``sequence_chunk``.
            NullSequenceException: If ``parent_or_seq_chunk_parent`` has no usable sequence ancestor.
            NoSuchAncestorException: If ``location`` has a ``sequence_chunk`` ancestor, but no ``chromosome`` ancestor.
                Such a relationship is required to lift from one chunk to a new chunk.
        """
        if parent_or_seq_chunk_parent is None:
            return location

        # if we are already a subset, we need to first lift back to genomic coordinates before lifting to this chunk
        if location.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            if not location.has_ancestor_of_type(SequenceType.CHROMOSOME):
                raise NoSuchAncestorException(
                    "This location does not have a chromosome ancestor of its sequence chunk, "
                    "which means it is not possible to lift to a new a chunk through the chromosome coordinates."
                )
            # ensure that both chromosomes are the same chromosome
            loc_chrom = location.first_ancestor_of_type(SequenceType.CHROMOSOME)
            par_chrom = parent_or_seq_chunk_parent.first_ancestor_of_type(SequenceType.CHROMOSOME)
            if loc_chrom.sequence and par_chrom.sequence:
                ObjectValidation.require_parents_equal_except_location(loc_chrom, par_chrom)
            else:
                ObjectValidation.require_parents_equal_except_location_and_sequence(loc_chrom, par_chrom)

            location = location.lift_over_to_first_ancestor_of_type(SequenceType.CHROMOSOME).reset_parent(
                parent_or_seq_chunk_parent.parent
            )

        if parent_or_seq_chunk_parent.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            if not parent_or_seq_chunk_parent.has_ancestor_of_type(SequenceType.CHROMOSOME):
                raise NoSuchAncestorException(
                    "Must have a chromosome in the hierarchy if a sequence chunk is provided."
                )

            chunk_parent = parent_or_seq_chunk_parent.first_ancestor_of_type(SequenceType.SEQUENCE_CHUNK)
            if not chunk_parent.sequence:
                raise NullSequenceException("Must have a sequence if a sequence chunk parent is provided.")

            location = location.reset_parent(chunk_parent.parent)
            sequence_chunk = chunk_parent.sequence
            # do not optimize blocks here -- this retains adjacent CDS intervals
            try:
                interval_location_rel_to_chunk = sequence_chunk.location_on_parent.parent_to_relative_location(
                    location, optimize_blocks=False
                )
            except LocationOverlapException:
                # the positions associated with this Location do not overlap the sequence chunk. However,
                # the chromosome location information can still be retained, but there is inherently no sequence
                # information.
                return EmptyLocation()
            interval_rel_to_chunk = interval_location_rel_to_chunk.reset_parent(parent_or_seq_chunk_parent)
            return interval_rel_to_chunk

        # since this is a whole genome (or something unknown), we don't need to lift anything up
        return location.reset_parent(parent_or_seq_chunk_parent)

    def liftover_to_parent_or_seq_chunk_parent(
        self,
        parent_or_seq_chunk_parent: Parent,
    ) -> "AbstractInterval":
        """
        This function returns a copy of this interval lifted over to a new coordinate system. If this interval
        is already in chunk-relative coordinates, it is first lifted back up the chromosome coordinates before
        the liftover occurs. This means that there *must* be a Parent somewhere in the ancestry with
        type "chromosome", and that Parent must match the supplied parent except for location information.

        Validation has to happen here in addition to in ``liftover_location_to_seq_chunk_parent()``, because
        at this point the parent of this current interval is still known. Once the ``to_dict()`` operation is performed,
        this information is list, and the new parent is applied under the assumption that it is valid.
        """
        if self.chunk_relative_location.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            if not self.chunk_relative_location.has_ancestor_of_type(SequenceType.CHROMOSOME):
                raise NoSuchAncestorException(
                    "This location does not have a chromosome ancestor of its sequence chunk, "
                    "which means it is not possible to lift to a new a chunk through the chromosome coordinates."
                )

        if self.chunk_relative_location.has_ancestor_of_type(SequenceType.CHROMOSOME):
            loc_chrom = self.chunk_relative_location.first_ancestor_of_type(SequenceType.CHROMOSOME)
            par_chrom = parent_or_seq_chunk_parent.first_ancestor_of_type(SequenceType.CHROMOSOME)
            if loc_chrom.sequence and par_chrom.sequence:
                ObjectValidation.require_parents_equal_except_location(loc_chrom, par_chrom)
            else:
                ObjectValidation.require_parents_equal_except_location_and_sequence(loc_chrom, par_chrom)

        return self.from_dict(self.to_dict(), parent_or_seq_chunk_parent)

    def _liftover_this_location_to_seq_chunk_parent(
        self,
        seq_chunk_parent: Parent,
    ):
        """Lift *this* interval to a new subset.

        This could happen as the result of a subsetting operation.

        This will introduce chunk-relative coordinates to this interval, or reduce the size of existing chunk-relative
        coordinates.

        This function calls the parent static method :meth:`AbstractInterval.liftover_location_to_seq_chunk_parent()`,
        but differs in two key ways:
        1. It acts on an instantiated subclass of this abstract class, modifying the location.
        2. It handles the case where a subclass is already a slice, by first lifting up to genomic coordinates.

        For these reasons, and particularly #1, this is a private method that is intended to be used during
        construction of a subclass. Modifying the locations in-place are generally a bad idea after initial
        construction of a interval class.
        """
        # if we are already a subset, we need to first lift back to genomic coordinates before lifting to this chunk
        if self._location.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            location = self._location.lift_over_to_first_ancestor_of_type(SequenceType.CHROMOSOME).reset_parent(
                seq_chunk_parent.parent
            )
        else:
            location = self._location
        self._location = self.liftover_location_to_seq_chunk_parent(location, seq_chunk_parent)

    def _reset_parent(self, parent: Optional[Parent] = None) -> None:
        """
        Convenience function that wraps location.reset_parent().

        NOTE: This function modifies this interval in-place, and does not return a new copy. This is different
        behavior than the base function, and is this way because this function is called recursively from collection
        objects.

        NOTE: Using this function presents the risk that you will change the sequence of this interval. There are no
        checks that the new parent provides the same sequence basis as the original parent.

        """
        self._location = self._location.reset_parent(parent)

    def has_ancestor_of_type(self, ancestor_type: Union[str, SequenceType]) -> bool:
        """
        Convenience function that wraps location.has_ancestor_of_type().
        """
        return self._location.has_ancestor_of_type(ancestor_type)

    def first_ancestor_of_type(self, ancestor_type: Union[str, SequenceType]) -> Parent:
        """
        Convenience function that returns the first ancestor of this type.
        """
        return self._location.first_ancestor_of_type(ancestor_type)

    def lift_over_to_first_ancestor_of_type(
        self, sequence_type: Optional[Union[str, SequenceType]] = SequenceType.CHROMOSOME
    ) -> Location:
        """
        Lifts the location member to another coordinate system. Is a no-op if there is no parent assigned.

        Returns:
            The lifted Location.
        """
        if self._location.parent is None:
            return self._location
        return self._location.lift_over_to_first_ancestor_of_type(sequence_type)

    def _import_qualifiers_from_list(self, qualifiers: Optional[Dict[Hashable, List[Hashable]]] = None):
        """Import input qualifiers to sets and store."""
        self.qualifiers = {}
        if qualifiers:
            if not isinstance(qualifiers, dict):
                raise ValidationException("Qualifiers must be a dictionary")
            for key, vals in qualifiers.items():
                if not isinstance(vals, list):
                    raise ValidationException("Qualifier values must be lists")
                self.qualifiers[key] = {str(x) for x in vals}

    def _export_qualifiers_to_list(self) -> Optional[Dict[Hashable, List[str]]]:
        """Export qualifiers back to lists. This is used when exporting to dictionary / converting back to marshmallow
        schemas.
        """
        if self.qualifiers:
            return {key: sorted(vals) for key, vals in self.qualifiers.items()}


class AbstractFeatureInterval(AbstractInterval, ABC):
    """This is a wrapper over :class:`~AbstractInterval` that adds functions shared across
    :class:`~biocantor.gene.transcript.TranscriptInterval`,
    :class:`~biocantor.gene.feature.FeatureInterval`,
    and :class:`~biocantor.gene.variants.VariantInterval`."""

    _genomic_ends: List[int]
    _genomic_starts: List[int]
    _strand: Strand
    _is_primary_feature: Optional[bool] = None

    def __len__(self):
        return sum((end - start) for end, start in zip(self._genomic_ends, self._genomic_starts))

    @lru_cache(maxsize=1)
    @property
    def chromosome_location(self) -> Location:
        """Returns the Location of this in *chromosome coordinates*.

        If the coordinate system is unknown, this will return the same coordinate system as
        ``chunk_relative_location``, that is the true underlying ``_location`` member.

        This Location object will always have the full span of the Interval in chromosome coordinates,
        even if this feature exists in chunk relative coordinates. As a result of this, if this
        Interval was built on chunk relative coordinates, the sequence information will not be present.
        """
        if self._parent_or_seq_chunk_parent and self._parent_or_seq_chunk_parent.has_ancestor_of_type(
            SequenceType.CHROMOSOME
        ):
            parent = self._parent_or_seq_chunk_parent.first_ancestor_of_type(SequenceType.CHROMOSOME)
            return CompoundInterval(self._genomic_starts, self._genomic_ends, self._strand, parent)
        else:
            return CompoundInterval(self._genomic_starts, self._genomic_ends, self._strand)

    @lru_cache(maxsize=1)
    @property
    def _chunk_relative_bounded_chromosome_location(self) -> Location:
        """
        Returns the Location of this in *chromosome coordinates*.

        This function is different from ``chromosome_location`` in that it will return a Location bounded
        by the chunk relative location of this Interval, if it exists.
        """
        if self.chunk_relative_location.is_empty:
            loc = CompoundInterval(self._genomic_starts, self._genomic_ends, self._strand)
            if self._parent_or_seq_chunk_parent.has_ancestor_of_type(SequenceType.CHROMOSOME):
                parent = self._parent_or_seq_chunk_parent.first_ancestor_of_type(SequenceType.CHROMOSOME)
                return loc.reset_parent(parent)
            else:
                return loc
        elif self.chunk_relative_location.has_ancestor_of_type(SequenceType.CHROMOSOME):
            return self.lift_over_to_first_ancestor_of_type(SequenceType.CHROMOSOME)
        return self.chunk_relative_location

    @lru_cache(maxsize=1)
    @property
    def chromosome_span(self) -> Location:
        """
        Returns the full span of this Interval in chromosome coordinates.
        """
        return self.chromosome_location._full_span_interval

    @lru_cache(maxsize=1)
    @property
    def chromosome_gaps_location(self) -> Location:
        """
        Returns the Location of the *gaps* of this Interval in chromosome coordinates. This is analogous to returning
        the intron coordinates.
        """
        return self.chromosome_location.gaps_location()

    @lru_cache(maxsize=1)
    @property
    def chunk_relative_span(self) -> Location:
        """
        Returns the full span of this Interval in chunk-relative coordinates.
        """
        return self.chunk_relative_location._full_span_interval

    @lru_cache(maxsize=1)
    @property
    def chunk_relative_gaps_location(self) -> Location:
        """
        Returns the Location of the *gaps* of this Interval in chunk-relative coordinates.
        This is analogous to returning the intron coordinates.
        """
        return self.chunk_relative_location.gaps_location()

    @abstractmethod
    def export_qualifiers(
        self, parent_qualifiers: Optional[Dict[Hashable, Set[Hashable]]] = None
    ) -> Dict[Hashable, Set[str]]:
        """Exports qualifiers for GFF3 or GenBank export. This merges top level keys with the arbitrary values"""

    @property
    def is_primary_feature(self) -> bool:
        """Is this the primary feature?"""
        return self._is_primary_feature is True

    @property
    def blocks(self) -> Iterable[SingleInterval]:
        """Wrapper for blocks function that reports blocks in chromosome coordinates"""
        yield from self.chromosome_location.blocks

    @property
    def relative_blocks(self) -> Iterable[SingleInterval]:
        """Wrapper for blocks function that reports blocks in chunk-relative coordinates"""
        yield from self.chunk_relative_location.blocks

    @abstractmethod
    def to_bed12(
        self,
        score: Optional[int] = 0,
        rgb: Optional[RGB] = RGB(0, 0, 0),
        name: Optional[str] = "feature_name",
        chromosome_relative_coordinates: bool = True,
    ) -> BED12:
        """Write a BED12 format representation of this :class:`AbstractFeatureInterval`.

        Both of these optional arguments are specific to the BED12 format.

        Args:
            score: An optional score associated with a interval. UCSC requires an integer between 0 and 1000.
            rgb: An optional RGB string for visualization on a browser. This allows you to have multiple colors
                on a single UCSC track.
            name: Which identifier in this record to use as 'name'. feature_name to guid. If the supplied string
                is not a valid attribute, it is used directly.
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.

        Return:
            A :class:`~biocantor.io.bed.BED12` object.

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
        """

    @abstractmethod
    def to_gff(
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict] = None,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterator[GFFRow]:
        """Writes a GFF format list of lists for this feature.

        The additional qualifiers are used when writing a hierarchical relationship back to files. GFF files
        are easier to work with if the children features have the qualifiers of their parents.

        Args:
            parent: ID of the Parent of this transcript.
            parent_qualifiers: Directly pull qualifiers in from this dictionary.
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
        """

    def sequence_pos_to_feature(self, pos: int) -> int:
        """Converts sequence position to relative position along this feature."""
        return self.chromosome_location.parent_to_relative_pos(pos)

    def sequence_interval_to_feature(self, chr_start: int, chr_end: int, chr_strand: Strand) -> Location:
        """Converts a contiguous interval on the sequence to a relative location within this feature."""
        loc = self.chromosome_location
        i = SingleInterval(chr_start, chr_end, chr_strand, parent=loc.parent)
        return loc.parent_to_relative_location(i)

    def feature_pos_to_sequence(self, pos: int) -> int:
        """Converts a relative position along this feature to sequence coordinate."""
        return self.chromosome_location.relative_to_parent_pos(pos)

    def feature_interval_to_sequence(self, rel_start: int, rel_end: int, rel_strand: Strand) -> Location:
        """Converts a contiguous interval relative to this feature to a spliced location on the sequence."""
        return self.chromosome_location.relative_interval_to_parent_location(rel_start, rel_end, rel_strand)

    def chunk_relative_pos_to_feature(self, pos: int) -> int:
        """Converts chunk-relative sequence position to relative position along this feature."""
        return self.chunk_relative_location.parent_to_relative_pos(pos)

    def chunk_relative_interval_to_feature(self, chr_start: int, chr_end: int, chr_strand: Strand) -> Location:
        """Converts a contiguous chunk-relative interval on the sequence to a relative location within this feature."""
        return self._location.parent_to_relative_location(
            SingleInterval(chr_start, chr_end, chr_strand, parent=self.chunk_relative_location.parent)
        )

    def feature_pos_to_chunk_relative(self, pos: int) -> int:
        """Converts a relative position along this feature to chunk-relative sequence coordinate."""
        return self.chunk_relative_location.relative_to_parent_pos(pos)

    def feature_interval_to_chunk_relative(self, rel_start: int, rel_end: int, rel_strand: Strand) -> Location:
        """
        Converts a contiguous interval relative to this feature to a chunk-relative spliced location on the sequence.
        """
        return self.chunk_relative_location.relative_interval_to_parent_location(rel_start, rel_end, rel_strand)

    @lru_cache(maxsize=1)
    def get_spliced_sequence(self) -> Sequence:
        """Returns the feature's *spliced*, *stranded* sequence."""
        ObjectValidation.require_location_has_parent_with_sequence(self._location)
        return self.chunk_relative_location.extract_sequence()

    @lru_cache(maxsize=1)
    def get_reference_sequence(self) -> Sequence:
        """Returns the feature's *unspliced*, *positive strand* genomic sequence."""
        ObjectValidation.require_location_has_parent_with_sequence(self._location)
        return self.chunk_relative_location.parent.sequence[self._location.start : self._location.end]

    @lru_cache(maxsize=1)
    def get_genomic_sequence(self) -> Sequence:
        """Returns the feature's *unspliced*, *stranded* (transcription orientation) genomic sequence."""
        ObjectValidation.require_location_has_parent_with_sequence(self._location)
        seq = self.chunk_relative_location.parent.sequence[self._location.start : self._location.end]
        if self.strand == Strand.PLUS:
            return seq
        else:
            return seq.reverse_complement()

    def _merge_qualifiers(
        self, other_qualifiers: Optional[Dict[Hashable, Set[str]]] = None
    ) -> Dict[Hashable, Set[str]]:
        """Merges this Interval's qualifiers dictionary with a new one, removing redundancy."""
        merged = self.qualifiers.copy()
        if other_qualifiers:
            for key, vals in other_qualifiers.items():
                if key not in merged:
                    merged[key] = set()
                merged[key].update(vals)
        return merged


class AbstractFeatureIntervalCollection(AbstractInterval, ABC):
    """
    Abstract class for holding groups of feature intervals. The two implementations of this class
    model Genes or non-transcribed FeatureCollections.

    These are always on the same sequence, but can be on different strands.
    """

    def __iter__(self):
        """Iterate over all children of this collection"""
        yield from self.iter_children()

    @property
    def strand(self) -> Strand:
        return self.chromosome_location.strand

    @abstractmethod
    def iter_children(self) -> Iterable["AbstractInterval"]:
        """Iterate over the children"""

    @abstractmethod
    def children_guids(self) -> Set[UUID]:
        """Get all of the GUIDs for children.

        Returns: A set of UUIDs
        """

    @abstractmethod
    def query_by_guids(self, id_or_ids: Union[UUID, List[UUID]]) -> "AbstractFeatureIntervalCollection":
        """Filter this collection object by a list of unique IDs.

        Args:
            id_or_ids: List of GUIDs, or unique IDs. Can also be a single ID.
        """

    def _reset_parent(self, parent: Optional[Parent] = None) -> None:
        """Reset parent of this collection, and all of its children.

        THIS FUNCTION IS ONLY INTENDED TO BE USED DURING INITIALIZATION OF A NEW INTERVAL OBJECT.
        USING THIS FUNCTION AFTER THAT POINT RUNS THE RISK OF THE PARENT OF THE OBJECT NOT BEING REFLECTED
        BY METHODS ON THIS FUNCTION THAT USE RESULT CACHING!

        NOTE: This function modifies this collection in-place, and does not return a new copy. This is different
        behavior than the base function, and is this way because all of the children of this collection are also
        recursively modified.

        NOTE: Using this function presents the risk that you will change the sequence of this interval. There are no
        checks that the new parent provides the same sequence basis as the original parent.

        This overrides :meth:`~biocantor.gene.feature.AbstractInterval.reset_parent()`. The original function
        will remain applied on the leaf nodes.
        """
        self._location = self._location.reset_parent(parent)
        for child in self:
            child._reset_parent(parent)

    def _initialize_location(self, start: int, end: int, parent_or_seq_chunk_parent: Optional[Parent] = None):
        """
        Initialize the location for this collection. Assumes that the start/end coordinates are genome-relative,
        and builds a chunk-relative location for this.

        Args:
            start: genome-relative start
            end: genome-relative end
            parent_or_seq_chunk_parent: A parent that could be null, genome relative, or sequence chunk relative.
        """
        self._location = SingleInterval(start, end, Strand.PLUS)
        if parent_or_seq_chunk_parent:
            if parent_or_seq_chunk_parent.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
                super()._liftover_this_location_to_seq_chunk_parent(parent_or_seq_chunk_parent)
            else:
                self._reset_parent(parent_or_seq_chunk_parent)

    def get_reference_sequence(self) -> Sequence:
        """Returns the *plus strand* sequence for this interval"""
        return self.chunk_relative_location.extract_sequence()

    @staticmethod
    def _find_primary_feature(
        intervals: Union[List["TranscriptInterval"], List["FeatureInterval"]]
    ) -> Optional[Union["TranscriptInterval", "FeatureInterval"]]:
        """
        Used in object construction to find the primary feature. Shared between :class:`GeneInterval`
        and :class:`FeatureIntervalCollection`.

        If not specified by the data source, primary features are determined by:

        1. If the feature is coding, then its CDS size
        2. The (spliced) feature size.
        3. The *position* of the feature within the ordered list of features.
        """
        # see if we were given a primary feature
        primary_feature = None
        for i, interval in enumerate(intervals):
            if interval.is_primary_feature:
                if primary_feature:
                    raise ValidationException("Multiple primary features/transcripts found")
                primary_feature = intervals[i]
        # if no primary interval was given, then infer by longest CDS then longest interval
        # if this is a feature, then there is no CDS, so set that value to 0
        if primary_feature is None:
            interval_sizes = sorted(
                (
                    [interval.cds_size if interval.interval_type == IntervalType.TRANSCRIPT else 0, len(interval), i]
                    for i, interval in enumerate(intervals)
                ),
                key=lambda x: (-x[0], -x[1]),
            )
            primary_feature = intervals[interval_sizes[0][2]]
        return primary_feature

    def _liftover_this_location_to_seq_chunk_parent(
        self,
        parent_or_seq_chunk_parent: Parent,
    ):
        """Lift over this collection and all of its children"""
        super()._liftover_this_location_to_seq_chunk_parent(parent_or_seq_chunk_parent)
        for child in self:
            child._liftover_this_location_to_seq_chunk_parent(parent_or_seq_chunk_parent)
