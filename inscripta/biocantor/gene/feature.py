"""
Object representation of features. Includes an abstract feature class that is also used by transcripts.

Each object is capable of exporting itself to BED and GFF3.
"""
from abc import ABC, abstractmethod
from typing import Optional, Any, Union, Dict, List, Set, Iterable, Hashable, TypeVar
from uuid import UUID

from inscripta.biocantor.exc import (
    EmptyLocationException,
    NoSuchAncestorException,
    ValidationException,
    NullSequenceException,
)
from inscripta.biocantor.gene.cds import CDSPhase
from inscripta.biocantor.io.bed import BED12, RGB
from inscripta.biocantor.io.gff3.constants import GFF_SOURCE, NULL_COLUMN, BioCantorFeatureTypes, BioCantorQualifiers
from inscripta.biocantor.io.gff3.exc import GFF3MissingSequenceNameError
from inscripta.biocantor.io.gff3.rows import GFFAttributes, GFFRow
from inscripta.biocantor.location.location import Location
from inscripta.biocantor.location.location_impl import SingleInterval, CompoundInterval
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent.parent import Parent
from inscripta.biocantor.sequence.sequence import Sequence
from inscripta.biocantor.util.bins import bins
from inscripta.biocantor.util.hashing import digest_object
from inscripta.biocantor.util.object_validation import ObjectValidation
from methodtools import lru_cache

# primitive data types possible as values of the list in a qualifiers dictionary
QualifierValue = TypeVar("QualifierValue", str, int, bool, float)


class AbstractInterval(ABC):
    """This is a wrapper over :class:`~biocantor.location.Location` that adds metadata coordinate transformation
    QOL functions."""

    location: Location
    _identifiers: List[Union[str, UUID]]
    qualifiers: Dict[Hashable, Set[str]]  # all subclasses convert qualifier values to sets of strings
    guid: UUID
    sequence_guid: Optional[UUID] = None
    sequence_name: Optional[str] = None
    bin: int

    def __len__(self):
        return len(self.location)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.to_dict() == other.to_dict()

    def __hash__(self):
        """Produces a hash, which is the GUID."""
        return hash(self.guid)

    @abstractmethod
    def to_dict(self, chromosome_relative_coordinates: bool = True) -> Dict[str, Any]:
        """Dictionary to build Model representation. Defaults to always exporting in original chromosome
        relative coordinates, but this can be disabled to export in sequence-chunk relative coordinates."""

    @staticmethod
    @abstractmethod
    def from_dict(vals: Dict[str, Any], parent_or_seq_chunk_parent: Optional[Parent] = None) -> "AbstractInterval":
        """Build an interval from a dictionary representation"""

    @abstractmethod
    def to_gff(
        self,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterable[GFFRow]:
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
    @abstractmethod
    def id(self) -> str:
        """Returns the ID of this feature. Provides a shared API across genes/transcripts and features."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Returns the name of this feature. Provides a shared API across genes/transcripts and features."""

    @property
    def start(self) -> int:
        """Returns genome relative start position."""
        return self.lift_over_to_first_ancestor_of_type("chromosome").start

    @property
    def end(self) -> int:
        """Returns genome relative end position."""
        return self.lift_over_to_first_ancestor_of_type("chromosome").end

    @property
    def chunk_relative_start(self) -> int:
        """Returns chunk relative start position."""
        return self.location.start

    @property
    def chunk_relative_end(self) -> int:
        """Returns chunk relative end position."""
        return self.location.end

    @property
    def strand(self) -> Strand:
        """Returns strand of location."""
        return self.location.strand

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
                    type="sequence_chunk",
                    parent=Parent(
                        location=SingleInterval(1, 15, Strand.PLUS,
                                               parent=Parent(id="genome_1_15", sequence_type="chromosome"))
                    ),
                )
            )

        Alternatively, if the sequence is coming straight from a file, it will be a :class:`Parent` with a
        :class:`Sequence` attached:

        .. code:: python

            parent = Parent(id="chr1", sequence=Sequence(genome, Alphabet.NT_STRICT, type="chromosome"))

        This convenience function detects which kind of parent is given, and sets up the appropriate location.

        Args:
            location: A location object, likely produced by :meth:`initialize_location()`.
            parent_or_seq_chunk_parent: An optional parent, either as a full chromosome or as a sequence chunk.

        Returns:
            A :class:`Location` object.

        Raises:
            ValidationException: If ``parent_or_seq_chunk_parent`` has no ancestor of type ``chromosome`` or
                ``sequence_chunk``.
            NullSequenceException: If ``parent_or_seq_chunk_parent`` has no usable sequence ancestor.
        """
        if parent_or_seq_chunk_parent is None:
            return location

        elif parent_or_seq_chunk_parent.has_ancestor_of_type("sequence_chunk"):

            chunk_parent = parent_or_seq_chunk_parent.first_ancestor_of_type("sequence_chunk")
            if not chunk_parent.sequence:
                raise NullSequenceException("Must have a sequence if a parent is provided.")

            location = location.reset_parent(chunk_parent.parent)
            sequence_chunk = chunk_parent.sequence
            interval_location_rel_to_chunk = sequence_chunk.location_on_parent.parent_to_relative_location(location)
            interval_rel_to_chunk = interval_location_rel_to_chunk.reset_parent(parent_or_seq_chunk_parent)

            return interval_rel_to_chunk

        # since this is a whole genome, we don't need to lift anything up
        elif parent_or_seq_chunk_parent.has_ancestor_of_type("chromosome"):
            if not parent_or_seq_chunk_parent.first_ancestor_of_type("chromosome").sequence:
                raise NullSequenceException("Must have a sequence if a parent is provided.")

            return location.reset_parent(parent_or_seq_chunk_parent)

        else:
            raise ValidationException("Provided Parent has no sequence of type 'chromosome' or 'sequence_chunk'.")

    def liftover_location_to_seq_chunk(
        self,
        seq_chunk_parent: Parent,
    ):
        """Lift this interval to a new subset.

        This could happen as the result of a subsetting operation.

        This will introduce chunk-relative coordinates to this interval, or reduce the size of existing chunk-relative
        coordinates.
        """
        if not seq_chunk_parent.has_ancestor_of_type("chromosome"):
            raise ValidationException("Provided Parent has no sequence of type 'chromosome'.")

        # if we are already a seq chunk, we need to lift ourselves back to genomic coordinates first
        if self.has_ancestor_of_type("sequence_chunk"):
            location = self.location.lift_over_to_first_ancestor_of_type("chromosome").reset_parent(
                seq_chunk_parent.parent
            )
        else:
            location = self.location.reset_parent(seq_chunk_parent.parent)
        sequence_chunk = seq_chunk_parent.sequence
        interval_location_rel_to_chunk = sequence_chunk.location_on_parent.parent_to_relative_location(location)
        interval_rel_to_chunk = interval_location_rel_to_chunk.reset_parent(seq_chunk_parent)

        self.location = interval_rel_to_chunk

    def reset_parent(self, parent: Parent):
        """
        Convenience function that wraps location.reset_parent().
        """
        self.location = self.location.reset_parent(parent)

    def has_ancestor_of_type(self, ancestor_type: str) -> bool:
        """
        Convenience function that wraps location.has_ancestor_of_type().
        """
        return self.location.has_ancestor_of_type(ancestor_type)

    def first_ancestor_of_type(self, ancestor_type: str) -> Parent:
        """
        Convenience function that returns the first ancestor of this type.
        """
        return self.location.first_ancestor_of_type(ancestor_type)

    def lift_over_to_first_ancestor_of_type(self, sequence_type: Optional[str] = "chromosome") -> Location:
        """
        Lifts the location member to another coordinate system. Is a no-op if there is no parent assigned.

        Returns:
            The lifted Location.
        """
        if self.location.parent is None:
            return self.location
        return self.location.lift_over_to_first_ancestor_of_type(sequence_type)

    def _import_qualifiers_from_list(self, qualifiers: Optional[Dict[Hashable, List[Hashable]]] = None):
        """Import input qualifiers to sets and store."""
        if qualifiers:
            self.qualifiers = {key: {str(x) for x in vals} for key, vals in qualifiers.items()}
        else:
            self.qualifiers = {}

    def _export_qualifiers_to_list(self) -> Optional[Dict[Hashable, List[str]]]:
        """Export qualifiers back to lists. This is used when exporting to dictionary / converting back to marshmallow
        schemas.
        """
        if self.qualifiers:
            return {key: sorted(vals) for key, vals in self.qualifiers.items()}


class AbstractFeatureInterval(AbstractInterval, ABC):
    """This is a wrapper over :class:`~AbstractInterval` that adds functions shared across
    :class:`~biocantor.gene.transcript.TranscriptInterval` and :class:`FeatureInterval`."""

    _is_primary_feature: Optional[bool] = None

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
        yield from self.lift_over_to_first_ancestor_of_type("chromosome").blocks

    @property
    def relative_blocks(self) -> Iterable[SingleInterval]:
        """Wrapper for blocks function that reports blocks in chunk-relative coordinates"""
        yield from self.location.blocks

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
    ) -> Iterable[GFFRow]:
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
        return self.lift_over_to_first_ancestor_of_type("chromosome").parent_to_relative_pos(pos)

    def sequence_interval_to_feature(self, chr_start: int, chr_end: int, chr_strand: Strand) -> Location:
        """Converts a contiguous interval on the sequence to a relative location within this feature."""
        loc = self.lift_over_to_first_ancestor_of_type("chromosome")
        i = SingleInterval(chr_start, chr_end, chr_strand, parent=loc.parent)
        return loc.parent_to_relative_location(i)

    def feature_pos_to_sequence(self, pos: int) -> int:
        """Converts a relative position along this feature to sequence coordinate."""
        return self.lift_over_to_first_ancestor_of_type("chromosome").relative_to_parent_pos(pos)

    def feature_interval_to_sequence(self, rel_start: int, rel_end: int, rel_strand: Strand) -> Location:
        """Converts a contiguous interval relative to this feature to a spliced location on the sequence."""
        return self.lift_over_to_first_ancestor_of_type("chromosome").relative_interval_to_parent_location(
            rel_start, rel_end, rel_strand
        )

    def chunk_relative_sequence_pos_to_feature(self, pos: int) -> int:
        """Converts chunk-relative sequence position to relative position along this feature."""
        return self.location.parent_to_relative_pos(pos)

    def chunk_relative_sequence_interval_to_feature(self, chr_start: int, chr_end: int, chr_strand: Strand) -> Location:
        """Converts a contiguous chunk-relative interval on the sequence to a relative location within this feature."""
        return self.location.parent_to_relative_location(
            SingleInterval(chr_start, chr_end, chr_strand, parent=self.location.parent)
        )

    def feature_pos_to_chunk_relative_sequence(self, pos: int) -> int:
        """Converts a relative position along this feature to chunk-relative sequence coordinate."""
        return self.location.relative_to_parent_pos(pos)

    def feature_interval_to_chunk_relative_sequence(self, rel_start: int, rel_end: int, rel_strand: Strand) -> Location:
        """
        Converts a contiguous interval relative to this feature to a chunk-relative spliced location on the sequence.
        """
        return self.location.relative_interval_to_parent_location(rel_start, rel_end, rel_strand)

    @lru_cache(maxsize=1)
    def get_spliced_sequence(self) -> Sequence:
        """Returns the feature's *spliced*, *stranded* sequence."""
        return self.location.extract_sequence()

    @lru_cache(maxsize=1)
    def get_reference_sequence(self) -> Sequence:
        """Returns the feature's *unspliced*, *positive strand* genomic sequence."""
        ObjectValidation.require_location_has_parent_with_sequence(self.location)
        return self.location.parent.sequence[self.location.start : self.location.end]

    @lru_cache(maxsize=1)
    def get_genomic_sequence(self) -> Sequence:
        """Returns the feature's *unspliced*, *stranded* (transcription orientation) genomic sequence."""
        seq = self.location.parent.sequence[self.location.start : self.location.end]
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


class FeatureInterval(AbstractFeatureInterval):
    """FeatureIntervals are generic intervals. These can be used to model genome promoters,
    open chromatin sites, etc.
    """

    _identifiers = ["feature_name", "feature_id"]

    def __init__(
        self,
        location: Location,
        qualifiers: Optional[Dict[Hashable, QualifierValue]] = None,
        sequence_guid: Optional[UUID] = None,
        sequence_name: Optional[str] = None,
        feature_types: Optional[List[str]] = None,
        feature_name: Optional[str] = None,
        feature_id: Optional[str] = None,
        guid: Optional[UUID] = None,
        feature_guid: Optional[UUID] = None,
        is_primary_feature: Optional[bool] = None,
    ):
        self.location = location
        self.sequence_guid = sequence_guid
        self.sequence_name = sequence_name
        self.feature_types = set(feature_types) if feature_types else set()  # stored as a set of types
        self.feature_name = feature_name
        self.feature_id = feature_id
        # qualifiers come in as a List, convert to Set
        self._import_qualifiers_from_list(qualifiers)
        self.bin = bins(self.start, self.end, fmt="bed")
        self._is_primary_feature = is_primary_feature

        if guid is None:
            self.guid = digest_object(
                self.location,
                self.qualifiers,
                self.sequence_name,
                self.feature_types,
                self.feature_name,
                self.feature_id,
                self.is_primary_feature,
            )
        else:
            self.guid = guid
        self.feature_guid = feature_guid

        if self.location.parent:
            ObjectValidation.require_location_has_parent_with_sequence(self.location)

    def __str__(self):
        return f"FeatureInterval(({self.location}), name={self.feature_name})"

    def __repr__(self):
        return "<{}>".format(str(self))

    @property
    def id(self) -> str:
        """Returns the ID of this feature. Provides a shared API across genes/transcripts and features."""
        return self.feature_id

    @property
    def name(self) -> str:
        """Returns the name of this feature. Provides a shared API across genes/transcripts and features."""
        return self.feature_name

    def to_dict(self, chromosome_relative_coordinates: bool = True) -> Dict[str, Any]:
        """Convert to a dict usable by :class:`biocantor.io.models.FeatureIntervalModel`."""
        blocks = self.blocks if chromosome_relative_coordinates else self.relative_blocks
        interval_starts, interval_ends = list(zip(*([x.start, x.end] for x in blocks)))
        return dict(
            interval_starts=interval_starts,
            interval_ends=interval_ends,
            strand=self.strand.name,
            qualifiers=self._export_qualifiers_to_list(),
            feature_id=self.feature_id,
            feature_name=self.feature_name,
            feature_types=sorted(self.feature_types) if self.feature_types else None,
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            feature_interval_guid=self.guid,
            feature_guid=self.feature_guid,
            is_primary_feature=self._is_primary_feature,
        )

    @staticmethod
    def from_dict(vals: Dict[str, Any], parent_or_seq_chunk_parent: Optional[Parent] = None) -> "FeatureInterval":
        """Build a :class:`FeatureInterval` from a dictionary."""
        location = FeatureInterval.initialize_location(
            vals["interval_starts"], vals["interval_ends"], Strand[vals["strand"]], parent_or_seq_chunk_parent
        )
        return FeatureInterval(
            location,
            qualifiers=vals["qualifiers"],
            sequence_guid=vals["sequence_guid"],
            sequence_name=vals["sequence_name"],
            feature_types=vals["feature_types"],
            feature_name=vals["feature_name"],
            feature_id=vals["feature_id"],
            guid=vals["feature_interval_guid"],
            feature_guid=vals["feature_guid"],
            is_primary_feature=vals["is_primary_feature"],
        )

    def intersect(
        self,
        location: Location,
        new_guid: Optional[UUID] = None,
        new_qualifiers: Optional[dict] = None,
    ) -> "FeatureInterval":
        """Returns a new FeatureInterval representing the intersection of this FeatureInterval's location with the
        other location.

        Strand of the other location is ignored; returned FeatureInterval is on the same strand as this FeatureInterval.
        """
        if not new_qualifiers:
            new_qualifiers = self.qualifiers

        location_same_strand = location.reset_strand(self.location.strand)
        intersection = self.location.intersection(location_same_strand)

        if intersection.is_empty:
            raise EmptyLocationException("Can't intersect disjoint intervals")

        return FeatureInterval(location=intersection, guid=new_guid, qualifiers=new_qualifiers)

    def export_qualifiers(
        self, parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None
    ) -> Dict[Hashable, Set[str]]:
        """Exports qualifiers for GFF3/GenBank export"""
        qualifiers = self._merge_qualifiers(parent_qualifiers)
        for key, val in [
            [BioCantorQualifiers.FEATURE_SYMBOL.value, self.feature_name],
            [BioCantorQualifiers.FEATURE_ID.value, self.feature_id],
        ]:
            if not val:
                continue
            if key not in qualifiers:
                qualifiers[key] = set()
            qualifiers[key].add(val)
        if self.feature_types:
            qualifiers[BioCantorQualifiers.FEATURE_TYPE.value] = self.feature_types
        return qualifiers

    def to_gff(
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterable[GFFRow]:
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
            GFF3MissingSequenceNameError: If there are no sequence names associated with this feature.
        """

        if not self.sequence_name:
            raise GFF3MissingSequenceNameError("Must have sequence names to export to GFF3.")

        if not chromosome_relative_coordinates and not self.has_ancestor_of_type("sequence_chunk"):
            raise NoSuchAncestorException(
                "Cannot export GFF in relative coordinates without a sequence_chunk ancestor."
            )

        qualifiers = self.export_qualifiers(parent_qualifiers)

        feature_id = str(self.guid)

        attributes = GFFAttributes(id=feature_id, qualifiers=qualifiers, name=self.feature_name, parent=parent)

        # "transcript" (feature interval) feature
        row = GFFRow(
            self.sequence_name,
            GFF_SOURCE,
            BioCantorFeatureTypes.FEATURE_INTERVAL,
            (self.start if chromosome_relative_coordinates else self.chunk_relative_start) + 1,
            self.end if chromosome_relative_coordinates else self.chunk_relative_end,
            NULL_COLUMN,
            self.strand,
            CDSPhase.NONE,
            attributes,
        )
        yield row

        # start adding exon features
        # re-use qualifiers, updating ID each time
        blocks = self.blocks if chromosome_relative_coordinates else self.relative_blocks
        for i, block in enumerate(blocks, 1):

            attributes = GFFAttributes(
                id=f"feature-{feature_id}-{i}", qualifiers=qualifiers, name=self.feature_name, parent=feature_id
            )
            row = GFFRow(
                self.sequence_name,
                GFF_SOURCE,
                BioCantorFeatureTypes.FEATURE_INTERVAL_REGION,
                block.start + 1,
                block.end,
                NULL_COLUMN,
                self.strand,
                CDSPhase.NONE,
                attributes,
            )
            yield row

    def to_bed12(
        self,
        score: Optional[int] = 0,
        rgb: Optional[RGB] = RGB(0, 0, 0),
        name: Optional[str] = "feature_name",
        chromosome_relative_coordinates: bool = True,
    ) -> BED12:
        """Write a BED12 format representation of this :class:`FeatureInterval`.

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
        # since blocks are iterated over twice, must be turned into a list otherwise the iterator is exhausted
        blocks = list(self.blocks) if chromosome_relative_coordinates else list(self.relative_blocks)
        block_sizes = [b.end - b.start for b in blocks]
        block_starts = [b.start - self.start for b in blocks]

        if chromosome_relative_coordinates:
            start = self.start
            end = self.end
        else:
            start = self.chunk_relative_start
            end = self.chunk_relative_end

        return BED12(
            self.sequence_name,
            start,
            end,
            getattr(self, name, name),
            score,
            self.strand,
            0,  # thickStart always 0 for non-coding
            0,  # thickEnd always 0 for non-coding
            rgb,
            self.location.num_blocks,
            block_sizes,
            block_starts,
        )
