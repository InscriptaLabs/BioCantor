import warnings
from itertools import count, zip_longest
from typing import Iterator, List, Union, Optional, Dict, Hashable, Any, Set, Tuple, TYPE_CHECKING, Type
from uuid import UUID

from methodtools import lru_cache

from biocantor.exc import (
    InvalidCDSIntervalError,
    NoSuchAncestorException,
    LocationOverlapException,
    MismatchedFrameException,
    EmptyLocationException,
)
from biocantor.gene.cds_frame import CDSPhase, CDSFrame
from biocantor.gene.codon import Codon, TranslationTable
from biocantor.gene.interval import AbstractFeatureInterval, QualifierValue
from biocantor.io.bed import RGB, BED12
from biocantor.io.gff3.constants import GFF_SOURCE, NULL_COLUMN, BioCantorFeatureTypes, BioCantorQualifiers
from biocantor.io.gff3.rows import GFFAttributes, GFFRow, GTFRow, GTFAttributes
from biocantor.location import Location, Strand, SingleInterval, CompoundInterval
from biocantor.parent import Parent, SequenceType
from biocantor.sequence import Sequence, Alphabet
from biocantor.util.hashing import digest_object

if TYPE_CHECKING:
    from biocantor.gene.variants import VariantIntervalCollection, VariantInterval


class CDSInterval(AbstractFeatureInterval):
    """
    This class represents a CDS interval, or an interval with coding potential. This is generally only used
    as a member of a :class:`~biocantor.gene.transcript.TranscriptInterval`. This class adds metadata and frame
    information to a Location object, and adds an understanding of codons, codon iteration, and translation.
    """

    frames = []
    _identifiers = ["protein_id", "product"]

    def __init__(
        self,
        cds_starts: List[int],
        cds_ends: List[int],
        strand: Strand,
        frames_or_phases: List[Union[CDSFrame, CDSPhase]],
        sequence_guid: Optional[UUID] = None,
        sequence_name: Optional[str] = None,
        protein_id: Optional[str] = None,
        product: Optional[str] = None,
        qualifiers: Optional[Dict[Hashable, QualifierValue]] = None,
        guid: Optional[UUID] = None,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ):
        self._location = self.initialize_location(cds_starts, cds_ends, strand, parent_or_seq_chunk_parent)
        self._genomic_starts = cds_starts
        self._genomic_ends = cds_ends
        self.start = cds_starts[0]
        self.end = cds_ends[-1]
        self._strand = strand
        self._parent_or_seq_chunk_parent = parent_or_seq_chunk_parent
        self.sequence_guid = sequence_guid
        self.sequence_name = sequence_name
        self.product = product
        self.protein_id = protein_id
        self._import_qualifiers_from_list(qualifiers)

        if len(frames_or_phases) != len(self._genomic_starts):
            raise MismatchedFrameException("Number of frame or phase entries must match number of exons")

        if len(self.chromosome_location) == 0:
            raise InvalidCDSIntervalError("Cannot have an empty CDS interval")

        # only allow either all CDSFrame or all CDSPhase
        is_frame = isinstance(frames_or_phases[0], CDSFrame)
        for frame_or_phase in frames_or_phases[1:]:
            if is_frame and isinstance(frame_or_phase, CDSPhase):
                raise MismatchedFrameException("Cannot mix frame and phase")
            elif not is_frame and isinstance(frame_or_phase, CDSFrame):
                raise MismatchedFrameException("Cannot mix frame and phase")

        if is_frame:
            self.frames = frames_or_phases
        else:
            self.frames = [x.to_frame() for x in frames_or_phases]

        if guid is None:
            self.guid = digest_object(
                self._genomic_starts,
                self._genomic_ends,
                self.strand,
                self.frames,
                self.product,
                self.protein_id,
                self.qualifiers,
            )
        else:
            self.guid = guid

        self._chunk_relative_codon_locations_cached = False

    def __str__(self):
        frame_str = ", ".join([str(p) for p in self.frames])
        return f"CDS(({self.chromosome_location}), ({frame_str})"

    def __repr__(self):
        return "<{}>".format(str(self))

    def __len__(self):
        return sum((end - start) for end, start in zip(self._genomic_ends, self._genomic_starts))

    @property
    def id(self) -> str:
        return self.protein_id

    @property
    def name(self) -> str:
        return self.product

    @property
    def chunk_relative_frames(self) -> List[CDSFrame]:
        """
        It may be the case that the chunk relative location of this CDSInterval object is a subset
        of the full chromosomal location. In this case, the frames list needs to be appropriately
        subsetted to the correct set of frame entries.

        However, it is far from trivial to subset frames in chunk context, because frames are calculated
        based on the full transcript length. Therefore, this function makes a blanket
        assumption that everything is in-frame within the interval of the chunk. In other words, if you are modeling
        a programmed frameshift using the Frames vector, this information will be lost. It does this by looping
        over the frames in transcription orientation until it finds the first exon that is within the chunk,
        then uses that to parameterize the frame generating function
        :meth:`CDSInterval.construct_frames_from_location()`
        """
        if not self.is_chunk_relative:
            return self.frames

        # TODO: Will there ever be a CDSFrame.NONE value? I don't think so because that is a genePred concept
        #   GFF3 will put Phase values only on CDS features, and GenBank is a lost cause
        fivep_phase = next(self._frame_iter(chunk_relative_frames=False)).to_phase().value
        distance_from_start = fivep_phase

        for genomic_exon in self._exon_iter(chunk_relative_exon=False):
            # chromosome location has overlapping blocks merged, so that the intersection always has one block
            # this is OK to do here since the original genomic intervals retain the overlapping information
            if isinstance(self._chunk_relative_bounded_chromosome_location, SingleInterval):
                chrom_loc = self._chunk_relative_bounded_chromosome_location
            else:
                chrom_loc = self._chunk_relative_bounded_chromosome_location.optimize_and_combine_blocks()

            # this is the section of the current exon with the chunk bounded location; it defines the
            # portion of this exon that is contained within the chunk
            intersection = genomic_exon.intersection(chrom_loc, match_strand=False)

            # this exon is entirely non-chunk
            if intersection.is_empty:
                distance_from_start += len(genomic_exon)
            # this should never happen, but just in case
            elif intersection.num_blocks != 1:
                raise LocationOverlapException("Found overlapping blocks after block optimization")
            # this exon is the first exon from the 5' end to be either fully or partially contained within the chunk
            # calculate the frames starting from here, after determining the distance from the 5' end
            else:
                if self.strand == Strand.PLUS:
                    distance_from_start += intersection.start - genomic_exon.start
                else:
                    distance_from_start += genomic_exon.end - intersection.end
                # this is equivalent to a CDSPhase, which we can convert to frame to get offset
                frame = CDSPhase(distance_from_start % 3).to_frame()
                return self.construct_frames_from_location(self.chunk_relative_location, frame)

    def to_dict(self, chromosome_relative_coordinates: bool = True) -> Dict[str, Any]:
        """
        Convert this CDS to a dictionary representation.

        If ``chromosome_relative_coordinates`` is ``False``, then the Frames list that comes out of this
        will lose programmed frameshift information.

        Args:
            chromosome_relative_coordinates: Optional flag to export the interval in chromosome relative
                or chunk-relative coordinates.

        Returns:
             A dictionary representation that can be passed to :meth:`CDSInterval.from_dict()`
        """
        if chromosome_relative_coordinates:
            cds_starts = self._genomic_starts
            cds_ends = self._genomic_ends
            cds_frames = [f.name for f in self.frames]
        else:
            cds_starts, cds_ends = list(zip(*([x.start, x.end] for x in self.chunk_relative_blocks)))
            cds_frames = [f.name for f in self.chunk_relative_frames]

        return dict(
            cds_starts=cds_starts,
            cds_ends=cds_ends,
            strand=self.strand.name,
            cds_frames=cds_frames,
            qualifiers=self._export_qualifiers_to_list(),
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            protein_id=self.protein_id,
            product=self.product,
        )

    @staticmethod
    def from_dict(vals: Dict[str, Any], parent_or_seq_chunk_parent: Optional[Parent] = None) -> "CDSInterval":
        """
        Construct a :class:`CDSInterval` from a dictionary representation such as one produced by
        :meth:`CDSInterval.to_dict()`.

        Args:
            vals: A dictionary representation.
            parent_or_seq_chunk_parent: An optional Parent to associate with this new interval.

        """
        return CDSInterval(
            cds_starts=vals["cds_starts"],
            cds_ends=vals["cds_ends"],
            strand=Strand[vals["strand"]],
            frames_or_phases=[CDSFrame[x] for x in vals["cds_frames"]],
            qualifiers=vals["qualifiers"],
            sequence_name=vals["sequence_name"],
            sequence_guid=vals["sequence_guid"],
            protein_id=vals["protein_id"],
            product=vals["product"],
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    @staticmethod
    def from_location(
        location: Location,
        cds_frames: List[Union[CDSFrame, CDSPhase]],
        sequence_guid: Optional[UUID] = None,
        sequence_name: Optional[str] = None,
        protein_id: Optional[str] = None,
        product: Optional[str] = None,
        qualifiers: Optional[Dict[Hashable, QualifierValue]] = None,
        guid: Optional[UUID] = None,
    ) -> "CDSInterval":
        """A convenience function that allows for construction of a :class:`CDSInterval` from a location object,
        a list of CDSFrames or CDSPhase, and optional metadata."""
        if location.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            raise NoSuchAncestorException(
                "Cannot call from_location with a chunk-relative location. Use from_chunk_relative_location()."
            )

        return CDSInterval(
            cds_starts=[x.start for x in location.blocks],
            cds_ends=[x.end for x in location.blocks],
            strand=location.strand,
            frames_or_phases=cds_frames,
            sequence_guid=sequence_guid,
            sequence_name=sequence_name,
            protein_id=protein_id,
            product=product,
            qualifiers=qualifiers,
            guid=guid,
            parent_or_seq_chunk_parent=location.parent,
        )

    @staticmethod
    def from_chunk_relative_location(
        location: Location,
        cds_frames: List[Union[CDSFrame, CDSPhase]],
        sequence_guid: Optional[UUID] = None,
        sequence_name: Optional[str] = None,
        protein_id: Optional[str] = None,
        product: Optional[str] = None,
        qualifiers: Optional[Dict[Hashable, QualifierValue]] = None,
        guid: Optional[UUID] = None,
    ) -> "CDSInterval":
        """
        Allows construction of a TranscriptInterval from a chunk-relative location. This is a location
        present on a sequence chunk, which should be built by the convenience function seq_chunk_to_parent:

        .. code-block:: python

            from biocantor.io.parser import seq_chunk_to_parent
            parent = seq_chunk_to_parent('AANAAATGGCGAGCACCTAACCCCCNCC', "NC_000913.3", 222213, 222241)
            loc = SingleInterval(5, 20, Strand.PLUS, parent=parent)

        And then, this can be lifted back to chromosomal coordinates like such:

        .. code-block:: python

            loc.lift_over_to_first_ancestor_of_type("chromosome")

        """
        if not location.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            raise NoSuchAncestorException("Must have a sequence chunk in the parent hierarchy.")

        chromosome_location = location.lift_over_to_first_ancestor_of_type(SequenceType.CHROMOSOME)
        return CDSInterval(
            cds_starts=[x.start for x in chromosome_location.blocks],
            cds_ends=[x.end for x in chromosome_location.blocks],
            strand=chromosome_location.strand,
            frames_or_phases=cds_frames,
            sequence_guid=sequence_guid,
            sequence_name=sequence_name,
            protein_id=protein_id,
            product=product,
            qualifiers=qualifiers,
            guid=guid,
            parent_or_seq_chunk_parent=location.parent,
        )

    def export_qualifiers(
        self, parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None
    ) -> Dict[Hashable, Set[Hashable]]:
        """Exports qualifiers for GFF3/GenBank export"""
        qualifiers = self._merge_qualifiers(parent_qualifiers)
        for key, val in [
            [BioCantorQualifiers.PROTEIN_ID.value, self.protein_id],
            [BioCantorQualifiers.PRODUCT.value, self.product],
        ]:
            if not val:
                continue
            if key not in qualifiers:
                qualifiers[key] = set()
            qualifiers[key].add(val)
        return qualifiers

    def _to_gff_or_gtf(
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None,
        chromosome_relative_coordinates: bool = True,
        raise_on_reserved_attributes: Optional[bool] = True,
        row_type: Union[Type[GFFRow], Type[GTFRow]] = GFFRow,
        attribute_type: Union[Type[GFFAttributes], Type[GTFAttributes]] = GFFAttributes,
    ) -> Iterator[Union[GFFRow, GTFRow]]:
        if not chromosome_relative_coordinates and not self.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            raise NoSuchAncestorException(
                "Cannot export GFF in relative coordinates without a sequence_chunk ancestor."
            )

        qualifiers = self.export_qualifiers(parent_qualifiers)

        cds_guid = str(self.guid)

        if chromosome_relative_coordinates:
            cds_blocks = zip(self._genomic_starts, self._genomic_ends)
            frames = self.frames
        else:
            cds_blocks = [[x.start, x.end] for x in self.chunk_relative_blocks]
            frames = self.chunk_relative_frames

        for i, block, frame in zip(count(1), cds_blocks, frames):
            start, end = block
            attributes = attribute_type(
                id=f"{cds_guid}-{i}",
                qualifiers=qualifiers,
                name=self.protein_id,
                parent=parent,
                raise_on_reserved_attributes=raise_on_reserved_attributes,
            )
            row = row_type(
                self.sequence_name,
                GFF_SOURCE,
                BioCantorFeatureTypes.CDS,
                start + 1,
                end,
                NULL_COLUMN,
                self.strand,
                frame.to_phase(),
                attributes,
            )
            yield row

    def to_gff(
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None,
        chromosome_relative_coordinates: bool = True,
        raise_on_reserved_attributes: Optional[bool] = True,
    ) -> Iterator[GFFRow]:
        """Writes a GFF format list of lists for this CDS.

        The additional qualifiers are used when writing a hierarchical relationship back to files. GFF files
        are easier to work with if the children features have the qualifiers of their parents.

        Args:
            parent: ID of the Parent of this transcript.
            parent_qualifiers: Directly pull qualifiers in from this dictionary.
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.
            raise_on_reserved_attributes: If ``True``, then GFF3 reserved attributes such as ``ID`` and ``Name`` present
                in the qualifiers will lead to an exception and not a warning.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
            GFF3MissingSequenceNameError: If there are no sequence names associated with this transcript.
        """
        yield from self._to_gff_or_gtf(
            parent,
            parent_qualifiers,
            chromosome_relative_coordinates,
            raise_on_reserved_attributes,
            GFFRow,
            GFFAttributes,
        )

    def to_gtf(
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterator[GTFRow]:
        """Writes a GTF format list of lists for this CDS.

        The additional qualifiers are used when writing a hierarchical relationship back to files. GTF files
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
            GFF3MissingSequenceNameError: If there are no sequence names associated with this transcript.
        """
        yield from self._to_gff_or_gtf(
            parent,
            parent_qualifiers,
            chromosome_relative_coordinates,
            False,
            GTFRow,
            GTFAttributes,
        )

    @property
    def has_canonical_start_codon(self) -> bool:
        """Does this CDS have a canonical valid start? Requires a sequence be associated."""
        return next(self.scan_codons()).is_canonical_start_codon

    def has_start_codon_in_specific_translation_table(
        self, translation_table: Optional[TranslationTable] = TranslationTable.DEFAULT
    ) -> bool:
        """
        Does this CDS have a valid start in a provided translation table? Requires a sequence be associated.

        Defaults to the ``DEFAULT`` table, which is just ``ATG``.
        """
        return next(self.scan_codons()).is_start_codon_in_specific_translation_table(translation_table)

    @property
    def has_valid_stop(self) -> bool:
        """Does this CDS have a valid stop? Requires a sequence be associated."""
        seq = self.extract_sequence()
        c = Codon(seq[-3:].sequence.upper())
        return c.is_stop_codon

    def _frame_iter(self, chunk_relative_frames: bool = True) -> Iterator[CDSFrame]:
        """Iterate over frames taking strand into account.

        If ``chunk_relative_frames`` is ``True``, then this iterator will only iterate over frames that
        overlap the relative chunk. These frames will potentially be reduced in quantity, and also shifted to handle
        exons that are now partial exons.
        """
        if chunk_relative_frames is True:
            if self.chunk_relative_location.strand == Strand.MINUS:
                yield from reversed(self.chunk_relative_frames)
            else:
                yield from self.chunk_relative_frames
        else:
            if self.strand == Strand.MINUS:
                yield from reversed(self.frames)
            else:
                yield from self.frames

    def _exon_iter(self, chunk_relative_exon: bool = True) -> Iterator[SingleInterval]:
        """Iterate over exons in transcription direction"""
        loc = self.chunk_relative_location if chunk_relative_exon else self.chromosome_location
        if loc.strand == Strand.PLUS or loc.strand == Strand.UNSTRANDED:
            yield from loc.blocks
        else:
            yield from reversed(loc.blocks)

    @lru_cache(maxsize=1)
    def extract_sequence(self) -> Sequence:
        """
        Returns a continuous CDS sequence that is in frame and always a multiple of 3.

        Any leading or trailing bases that are annotated as CDS but cannot form a full codon
        are removed. Additionally, any internal codons that are incomplete are removed.

        Incomplete internal codons are determined by comparing the CDSFrame of each exon
        as annotated, to the expected value of the CDSFrame. This allows for an annotation
        to model things like programmed frameshifts and indels that may be assembly errors.

        This function has been optimized to run as fast as possible. The original implementation iterated
        over every codon, but this is slower because it has a lot of object instantiation overhead. However,
        if those objects have already been instantiated and cached, then it is faster to just re-use them.
        """
        if self._chunk_relative_codon_locations_cached is True:
            codons = (str(codon_location.extract_sequence()) for codon_location in self.chunk_relative_codon_locations)
            seq = "".join(codons)
            return seq
        if self.num_blocks > 1:
            window_fn = self._prepare_multi_exon_window_for_scan_codon_locations
        else:
            window_fn = self._prepare_single_exon_window_for_scan_codon_locations
        location, offset = window_fn(relative_window=None, chunk_relative_coordinates=True)
        seq = str(location.extract_sequence())[offset : len(location) - ((len(location) - offset) % 3)]
        return Sequence(seq, Alphabet.NT_EXTENDED, validate_alphabet=False)

    @property
    def num_codons(self) -> int:
        """
        Returns the total number of codons. This will reflect the true number of codons,
        even if this CDSInterval is parented on a sequence chunk.
        """
        return len(self.chromosome_codon_locations)

    @property
    def num_chunk_relative_codons(self) -> int:
        """
        Returns the number of codons.

        NOTE: If this CDS is a subset of the original sequence, this number will represent the subset,
        not the original size!

        Any leading or trailing bases that are annotated as CDS but cannot form a full codon
        are excluded. Additionally, any internal codons that are incomplete are excluded.

        Incomplete internal codons are determined by comparing the CDSFrame of each exon
        as annotated, to the expected value of the CDSFrame. This allows for an annotation
        to model things like programmed frameshifts and indels that may be assembly errors.
        """
        return len(self.chunk_relative_codon_locations)

    def scan_codons(self, truncate_at_in_frame_stop: Optional[bool] = False) -> Iterator[Codon]:
        """
        Iterator along codons. If truncate_at_in_frame_stop is True,
        this will stop iteration at the first in-frame stop.
        """
        seq = self.extract_sequence()
        for i in range(0, len(seq), 3):
            c = Codon(str(seq[i : i + 3]).upper())
            yield c
            if truncate_at_in_frame_stop and c.is_stop_codon:
                break

    @lru_cache(maxsize=1)
    @property
    def chunk_relative_codon_locations(self) -> Tuple[Location]:
        """
        Returns a tuple of codon locations in *chunk relative* coordinates.

        This function calls ``scan_codon_locations`` and stores the full result as a cached
        tuple.
        """
        self._chunk_relative_codon_locations_cached = True
        return tuple(self.scan_chunk_relative_codon_locations())

    @lru_cache(maxsize=1)
    @property
    def chromosome_codon_locations(self) -> Tuple[Location]:
        """
        Returns a tuple of codon locations in *chromosome* coordinates.

        If this is a chunk-relative CDS, the returned locations will not have sequence information.

        This function calls ``scan_codon_locations`` and stores the full result as a cached
        tuple.
        """
        return tuple(self._scan_codon_locations(chunk_relative_coordinates=False))

    def _convert_chromosome_start_end_to_relative_window(
        self,
        chromosome_start: Optional[int] = None,
        chromosome_end: Optional[int] = None,
        expand_window_to_partial_codons: bool = False,
    ) -> Optional[SingleInterval]:
        """
        Converts possibly null chromosomal start/end positions to a SingleInterval representing the genomic
        span of that window. Null start/end values default to the start/end of this CDS.
        """
        if chromosome_start is None and chromosome_end is None:
            return None
        if chromosome_start is None:
            chromosome_start = self.start
        if chromosome_end is None:
            chromosome_end = self.end
        if expand_window_to_partial_codons:
            chromosome_start, chromosome_end = self._expand_coordinates_to_codons(chromosome_start, chromosome_end)
        relative_window = SingleInterval(chromosome_start, chromosome_end, self.strand, self.chromosome_location.parent)
        return relative_window

    def _expand_coordinates_to_codons(self, chromosome_start: int, chromosome_end: int) -> Tuple[int, int]:
        """
        Convenience function to take a pair of chromosome coordinates and return new coordinates that contain only
        full codons.
        """

        if chromosome_start < self.chromosome_location.start:
            chromosome_start = self.start
        if chromosome_end > self.chromosome_location.end:
            chromosome_end = self.end
        cds_interval = self.sequence_interval_to_cds(chromosome_start, chromosome_end, Strand.PLUS)
        adjusted_cds_start = cds_interval.start - (cds_interval.start % 3)
        adjusted_cds_end = cds_interval.end - (cds_interval.end % -3)
        chromosome_interval = self.cds_interval_to_sequence(adjusted_cds_start, adjusted_cds_end, Strand.PLUS)
        return chromosome_interval.start, chromosome_interval.end

    def scan_chunk_relative_codon_locations(
        self,
        chromosome_start: Optional[int] = None,
        chromosome_end: Optional[int] = None,
        expand_window_to_partial_codons: bool = False,
    ) -> Iterator[Location]:
        """
        Returns an iterator over codon locations in *chunk relative* coordinates.

        Any leading or trailing bases that are annotated as CDS but cannot form a full codon
        are excluded. Additionally, any internal codons that are incomplete are excluded.

        Incomplete internal codons are determined by comparing the CDSFrame of each exon
        as annotated, to the expected value of the CDSFrame. This allows for an annotation
        to model things like programmed frameshifts and indels that may be assembly errors.

        Args:
            chromosome_start: An optional *chromosome* position to offset the iteration to. The resulting codons
                will maintain frame.
            chromosome_end: An optional *chromosome* position to offset the iteration to end at. The resulting codons
                will maintain frame. This number can be larger than the chromosome end position.
            expand_window_to_partial_codons: If ``True``, and the ``chromosome_start`` or ``chromosome_end`` parameters
                are set to values within a codon, the full codon will be retained. If ``False``, partial codons
                will be eliminated.
        """
        yield from self._scan_codon_locations(
            self._convert_chromosome_start_end_to_relative_window(
                chromosome_start, chromosome_end, expand_window_to_partial_codons
            ),
            chunk_relative_coordinates=True,
        )

    def scan_chromosome_codon_locations(
        self,
        chromosome_start: Optional[int] = None,
        chromosome_end: Optional[int] = None,
        expand_window_to_partial_codons: bool = False,
    ) -> Iterator[Location]:
        """
        Returns an iterator over codon locations in *chromosome* coordinates.

        If this is a chunk-relative CDS, the returned locations will not have sequence information.

        Any leading or trailing bases that are annotated as CDS but cannot form a full codon
        are excluded. Additionally, any internal codons that are incomplete are excluded.

        Incomplete internal codons are determined by comparing the CDSFrame of each exon
        as annotated, to the expected value of the CDSFrame. This allows for an annotation
        to model things like programmed frameshifts and indels that may be assembly errors.

        Args:
            chromosome_start: An optional *chromosome* position to offset the iteration to. The resulting codons
                will maintain frame.
            chromosome_end: An optional *chromosome* position to offset the iteration to end at. The resulting codons
                will maintain frame. This number can be larger than the chromosome end position.
            expand_window_to_partial_codons: If ``True``, and the ``chromosome_start`` or ``chromosome_end`` parameters
                are set to values within a codon, the full codon will be retained. If ``False``, partial codons
                will be eliminated.
        """
        yield from self._scan_codon_locations(
            self._convert_chromosome_start_end_to_relative_window(
                chromosome_start, chromosome_end, expand_window_to_partial_codons
            ),
            chunk_relative_coordinates=False,
        )

    def scan_codon_locations(self) -> Iterator[Location]:
        """
        Scan codon locations in chunk-relative coordinates. This function exists for backwards compatibility
        and is deprecated. It however retains the speedup optimization introduced in BioCantor 0.10.0.
        """
        warnings.warn(
            DeprecationWarning(
                "CDSInterval.scan_codon_locations is deprecated. "
                "Use scan_chromosome_codon_locations or scan_chunk_relative_codon_locations."
            )
        )
        yield from self._scan_codon_locations(chunk_relative_coordinates=True)

    def _scan_codon_locations(
        self, relative_window: Optional[SingleInterval] = None, chunk_relative_coordinates: bool = True
    ) -> Iterator[Location]:
        """
        Returns an iterator over codon locations.

        ``relative_window`` must be in *chromosome* coordinates, ideally built by
        ``_convert_chromosome_start_end_to_relative_window``.

        Any leading or trailing bases that are annotated as CDS but cannot form a full codon
        are excluded.
        """
        # can only do naive window scanning if this CDS has exactly one exon
        if self.num_blocks > 1:
            codon_fn = self._prepare_multi_exon_window_for_scan_codon_locations
        else:
            codon_fn = self._prepare_single_exon_window_for_scan_codon_locations
        location, offset = codon_fn(relative_window, chunk_relative_coordinates)
        # must make sure this CDS (with its offset) is at least one codon long
        if (len(location) - offset) >= 3:
            yield from location.scan_windows(3, 3, offset)

    @lru_cache(maxsize=20)
    def _prepare_single_exon_window_for_scan_codon_locations(
        self, relative_window: Optional[SingleInterval] = None, chunk_relative_coordinates: bool = True
    ) -> Tuple[Location, int]:
        """
        This function exists to prepare a Location object to pass to the iterator ``_scan_codon_locations``.
        By placing the logic in this function, the result can be cached, as you cannot cache iterators.

        ``relative_window`` must be in *chromosome* coordinates, ideally built by
        ``_convert_chromosome_start_end_to_relative_window``.

        Returns a tuple of the Location to be iterated over, and its offset.
        """
        # do all initial work in chromosome coordinates
        loc = self.chromosome_location
        offset = self.frames[0].value

        if relative_window:
            relative_loc = loc.intersection(relative_window)
        else:
            relative_loc = loc

        if chunk_relative_coordinates and self.is_chunk_relative:
            # lift the cleaned window on to chunk relative coordinate system
            chunk_relative_cleaned_location = self.liftover_location_to_seq_chunk_parent(
                relative_loc, self.chunk_relative_location.parent
            )
            # lift this back to chromosome coordinates -- this produces a chromosome coordinate Location
            # whose bounds are the portion of this CDS that are contained on the sequence chunk
            loc_on_chrom = chunk_relative_cleaned_location.lift_over_to_first_ancestor_of_type(SequenceType.CHROMOSOME)
            offset += self._calculate_frame_offset(relative_loc, loc_on_chrom)
            return chunk_relative_cleaned_location, offset
        else:
            offset += self._calculate_frame_offset(loc, relative_loc)
            return relative_loc, offset

    @lru_cache(maxsize=20)
    def _prepare_multi_exon_window_for_scan_codon_locations(
        self, relative_window: Optional[SingleInterval] = None, chunk_relative_coordinates: bool = True
    ) -> Tuple[Location, int]:
        """
        This function exists to prepare a Location object to pass to the iterator ``_scan_codon_locations``.
        By placing the logic in this function, the result can be cached, as you cannot cache iterators.

        Returns a tuple of the Location to be iterated over, and its offset.
        """
        next_frame = CDSFrame.ZERO
        cleaned_rel_starts = []
        cleaned_rel_ends = []
        loc = self.chromosome_location
        # zip_longest is used here to ensure that the two iterators are always actually in sync
        for exon, frame in zip_longest(self._exon_iter(False), self._frame_iter(False)):
            if exon is None or frame is None:
                raise MismatchedFrameException("Frame iterator is not in sync with exon iterator")

            start_to_rel = loc.parent_to_relative_pos(exon.start)
            end_to_rel_inclusive = loc.parent_to_relative_pos(exon.end - 1)
            rel_start = min(start_to_rel, end_to_rel_inclusive)
            rel_end = max(start_to_rel, end_to_rel_inclusive) + 1
            if next_frame != frame:
                rel_start += frame.value
                # remove trailing codon from previous block
                shift = sum((coords[1] - coords[0] for coords in zip(cleaned_rel_starts, cleaned_rel_ends))) % 3
                if shift > 0:
                    # it may be possible for this shift to end up producing a 0bp block
                    # this will be dropped in the list comprehension below that generates the cleaned_blocks
                    cleaned_rel_ends[-1] = cleaned_rel_ends[-1] - shift
                # we are now inherently in frame
                next_frame = CDSFrame.ZERO
            # it may be the case that the removal of the trailing codon from the previous block entirely
            # eliminates that block. In that case, skip it.
            if rel_start >= rel_end:
                continue
            cleaned_rel_starts.append(rel_start)
            cleaned_rel_ends.append(rel_end)
            # this is what we expect the next frame to be, if no frameshift occurred
            next_frame = next_frame.shift(rel_end - rel_start)
        cleaned_blocks = [
            loc.relative_interval_to_parent_location(cleaned_rel_starts[i], cleaned_rel_ends[i], Strand.PLUS)
            for i in range(len(cleaned_rel_starts))
            # remove 0bp blocks here to avoid having to call optimize_blocks()
            if cleaned_rel_ends[i] != cleaned_rel_starts[i]
        ]
        cleaned_location = CompoundInterval.from_single_intervals(cleaned_blocks)

        if relative_window:
            relative_cleaned_location = cleaned_location.intersection(relative_window)
        else:
            relative_cleaned_location = cleaned_location

        if chunk_relative_coordinates and self.is_chunk_relative:
            # lift the cleaned window on to chunk relative coordinate system
            chunk_relative_cleaned_location = self.liftover_location_to_seq_chunk_parent(
                relative_cleaned_location, self.chunk_relative_location.parent
            )
            # lift this back to chromosome coordinates -- this produces a chromosome coordinate Location
            # whose bounds are the portion of this CDS that are contained on the sequence chunk
            loc_on_chrom = chunk_relative_cleaned_location.lift_over_to_first_ancestor_of_type(SequenceType.CHROMOSOME)
            offset = self._calculate_frame_offset(relative_cleaned_location, loc_on_chrom)
            return chunk_relative_cleaned_location, offset
        else:
            offset = self._calculate_frame_offset(cleaned_location, relative_cleaned_location)
            return relative_cleaned_location, offset

    def _calculate_frame_offset(self, cleaned_location: Location, loc_on_chrom: Location) -> int:
        """
        In either single-exon or multi-exon codon iteration, if this CDSInterval exists on chunk-relative coordinates
        that slice down the CDS, then the initial offset provided by the CDSFrame field must be adjusted
        to maintain frame.
        """
        # determine the offset needed for codon iteration to be in-frame
        # calculate the number of bases from the 5' end that this sequence chunk removes, if any
        if self.strand == Strand.PLUS:
            fivep_loc = cleaned_location.relative_interval_to_parent_location(
                0, cleaned_location.parent_to_relative_pos(loc_on_chrom.start), Strand.PLUS
            )
        else:
            fivep_loc = cleaned_location.relative_interval_to_parent_location(
                0, cleaned_location.parent_to_relative_pos(loc_on_chrom.end - 1), Strand.PLUS
            )
        # this is equivalent to a CDSPhase, which we can convert to frame to get offset
        fivep_distance_mod3 = len(fivep_loc) % 3
        phase = CDSPhase(fivep_distance_mod3)
        offset = phase.to_frame().value
        return offset

    @lru_cache(maxsize=2)
    def translate(
        self,
        truncate_at_in_frame_stop: Optional[bool] = False,
        translation_table: Optional[TranslationTable] = TranslationTable.DEFAULT,
        strict: bool = True,
    ) -> Sequence:
        """
        Returns amino acid sequence of this CDS.

        Parameters
        ----------
        truncate_at_in_frame_stop
            If truncate_at_in_frame_stop is ``True``, this will stop at the first in-frame stop.
        translation_table
            Currently the ``translation_table`` field only controls the start codon. Using non-standard
            translation tables will change the set of start codons that code for Methionine,
            and will not change any other codons.
        strict
            If False, allows untranslatable codons to be represented by an X. Otherwise, throws ValueError.

        Returns
        -------
        Sequence
            The translated amino acid sequence

        Raises
        ------
        ValueError
            Codon is untranslatable and allow_unknown_translation is False
        """
        seq = str(self.extract_sequence()).upper()
        translated_seq = []
        for i in range(0, len(seq), 3):
            codon_str = seq[i : i + 3]

            codon = Codon(codon_str)
            if i == 0 and codon.is_start_codon_in_specific_translation_table(translation_table):
                translated_seq.append(Codon("ATG").translate())
            else:
                if strict and not codon.is_strict_codon:
                    raise ValueError(f"Codon is not a strict codon: '{codon}'")
                translated_seq.append(codon.translate(strict=strict))

            if truncate_at_in_frame_stop and codon.is_stop_codon and i != len(seq) - 3:
                break
        alphabet = Alphabet.AA if strict else Alphabet.AA_STRICT_UNKNOWN
        return Sequence("".join(translated_seq), alphabet, validate_alphabet=False)

    @lru_cache(maxsize=1)
    @property
    def has_in_frame_stop(self) -> bool:
        """Does this CDS have an in-frame stop codon?"""
        return "*" in str(self.translate()[:-1])

    @staticmethod
    def construct_frames_from_location(
        location: Location, starting_frame: Optional[CDSFrame] = CDSFrame.ZERO
    ) -> List[CDSFrame]:
        """
            Construct a list of CDSFrames from a Location. This is intended to construct frames in situations where
            the frames are not known. One example of such a case is when parsing GenBank files, which have only
            a ``codon_start`` field to measure the offset at the start of translation.

            This function is extremely hard to understand, so I hope the below example helps:


            1. Plus strand:

            .. code-block::

                CompoundInterval([0, 7, 12], [5, 11, 18], Strand.PLUS)


            .. code-block::

                Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
                Sequence:   A A A C A A A A G G G  T  A  C  C  C  A  A  A  A  A  A
                Exons:      A A A C A     A G G G     A  C  C  C  A  A
                Zero Frame: 0 1 2 0 1     2 0 1 2     0  1  2  0  1  2
                One Frame:  - 0 1 2 0     1 2 0 1     2  0  1  2  0  1
                Two Frame:  - - 0 1 2     0 1 2 0     1  2  0  1  2  0


            In the non-zero case, the ``[0, 1, 2]`` cycle is offset by 1 or 2 bases.

            So, for this test case we expect the frames to be:

        .. code-block::

                Zero Frame: [0, 2, 0]
                One Frame:  [1, 1, 2]
                Two Frame:  [2, 0, 1]


            2. Minus strand:

            .. code-block::

                CompoundInterval([0, 7, 12], [5, 11, 18], Strand.MINUS)


            .. code-block::

                Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
                Sequence:   A A A C A A A A G G G  T  A  C  C  C  A  A  A  A  A  A
                Exons:      A A A C A     A G G G     A  C  C  C  A  A
                Zero Frame: 2 1 0 2 1     0 2 1 0     2  1  0  2  1  0
                One Frame:  1 0 2 1 0     2 1 0 2     1  0  2  1  0  -
                Two Frame:  0 2 1 0 2     1 0 2 1     0  2  1  0  -  -


            Now, for negative strand CDS intervals, the frame list is still in plus strand orientation.

            So, for this test case we expect the frames to be:

            .. code-block::

                Zero Frame: [1, 0, 0]
                One Frame:  [0, 2, 1]
                Two Frame:  [2, 1, 2]

            Args:
                location: A interval of the CDS.
                starting_frame: Frame to start iteration with. If ``codon_start`` was the source of this value,
                    then you would subtract one before converting to :class:`CDSFrame`.

            Returns:
                A list of :class:`CDSFrame` that could be combined with the input Location to build a
                :class:`CDSInterval`.
        """
        # edge case: if there is only one block, then just return the starting frame
        if location.num_blocks == 1:
            return [starting_frame]
        # find size of every block except last block in transcription orientation
        sizes = [len(x) for x in location.scan_blocks()][:-1]
        # shift first block size to starting frame
        sizes[0] -= starting_frame.value
        # start in frame with new shifted block size
        frames = [CDSFrame.ZERO]
        for s in sizes:
            frames.append(frames[-1].shift(s))
        # swap back in original starting frame
        frames[0] = starting_frame
        # flip around if this is negative strand because frames are in + orientation
        if location.strand == Strand.MINUS:
            frames = frames[::-1]
        return frames

    def optimize_blocks(self) -> "CDSInterval":
        """
        Combine the blocks of this CDS interval, preserving overlapping blocks.

        Once this operation is performed, internal frameshifts modeled by 0bp gaps will be lost, and the resulting
        translation will be out of frame downstream.

        Returns:
            A new :class:`CDSInterval` that has been merged.
        """
        new_loc = self.chunk_relative_location.optimize_blocks()
        first_frame = next(self._frame_iter())
        frames = CDSInterval.construct_frames_from_location(new_loc, first_frame)
        return CDSInterval.from_location(new_loc, frames)

    def optimize_and_combine_blocks(self) -> "CDSInterval":
        """
        Combine the blocks of this CDS interval, including removing overlapping blocks.

        Once this operation is performed, internal frameshifts modeled by 0bp gaps will be lost, as well as programmed
        frameshifts modeled by overlapping blocks. The resulting translations will be out of frame downstream.

        Returns:
            A new :class:`CDSInterval` that has been merged.
        """
        if isinstance(self.chunk_relative_location, CompoundInterval):
            new_loc = self.chunk_relative_location.optimize_and_combine_blocks()
        else:
            new_loc = self.chunk_relative_location
        first_frame = next(self._frame_iter())
        frames = CDSInterval.construct_frames_from_location(new_loc, first_frame)
        return CDSInterval.from_location(new_loc, frames)

    def to_bed12(
        self,
        score: Optional[int] = 0,
        rgb: Optional[RGB] = RGB(0, 0, 0),
        name: Optional[str] = "feature_name",
        chromosome_relative_coordinates: bool = True,
    ) -> BED12:
        raise NotImplementedError

    def cds_pos_to_sequence(self, pos: int) -> int:
        """Converts a relative position along the CDS to sequence coordinate."""
        return self.chromosome_location.relative_to_parent_pos(pos)

    def cds_pos_to_chunk_relative(self, pos: int) -> int:
        """Converts a relative position along the CDS to chunk-relative sequence coordinate."""
        return self.chunk_relative_location.relative_to_parent_pos(pos)

    def cds_interval_to_sequence(self, rel_start: int, rel_end: int, rel_strand: Strand) -> Location:
        """Converts a contiguous interval relative to the CDS to a spliced location on the sequence."""
        return self.chromosome_location.relative_interval_to_parent_location(rel_start, rel_end, rel_strand)

    def cds_interval_to_chunk_relative(self, rel_start: int, rel_end: int, rel_strand: Strand) -> Location:
        """Converts a contiguous interval relative to the CDS to a spliced location on the chunk-relative sequence."""
        return self.chunk_relative_location.relative_interval_to_parent_location(rel_start, rel_end, rel_strand)

    def sequence_pos_to_cds(self, pos: int) -> int:
        """
        Converts a *sequence* relative position to a CDS position. This is the distance from the translation
        start in CDS coordinates.

        Returns:
            An integer position in CDS coordinates.
        Raises:
            InvalidPositionException: If the position provided is not part of this CDSInterval.
        """
        return self.chromosome_location.parent_to_relative_pos(pos)

    def chunk_relative_pos_to_cds(self, pos: int) -> int:
        """Converts chunk-relative sequence position to relative position along the CDS."""
        return self.chunk_relative_location.parent_to_relative_pos(pos)

    def sequence_interval_to_cds(self, chr_start: int, chr_end: int, chr_strand: Strand) -> Location:
        """Converts a contiguous interval on the sequence to a relative location within the CDS."""
        i = SingleInterval(chr_start, chr_end, chr_strand, parent=self.chromosome_location.parent)
        return self.chromosome_location.parent_to_relative_location(i)

    def chunk_relative_interval_to_cds(self, chr_start: int, chr_end: int, chr_strand: Strand) -> Location:
        """Converts a contiguous interval on the chunk-relative sequence to a relative location within the CDS."""
        return self.chunk_relative_location.parent_to_relative_location(
            SingleInterval(chr_start, chr_end, chr_strand, parent=self.chunk_relative_location.parent)
        )

    def sequence_pos_to_amino_acid(self, pos: int) -> int:
        """
        Converts a *sequence* relative position to amino acid. The resulting value is always left-aligned.

        Returns:
            A zero based integer position on the amino acid sequence, left aligned.
        Raises:
            InvalidPositionException: If the position provided is not part of this CDSInterval.
        """
        return self.sequence_pos_to_cds(pos) // 3

    def incorporate_variants(self, variants: Union["VariantInterval", "VariantIntervalCollection"]) -> "CDSInterval":
        """
        Incorporate all of the variant(s) for an input VariantInterval or VariantIntervalCollection,
        producing a new CDSInterval with those changes incorporated.
        """
        new_loc = variants.lift_over_location(self.chunk_relative_location)
        if new_loc.is_empty:
            raise EmptyLocationException("Variant incorporation led to an EmptyLocation")
        fn = CDSInterval.from_chunk_relative_location if self.is_chunk_relative else CDSInterval.from_location
        new_frames = CDSInterval.construct_frames_from_location(new_loc, self.frames[0])
        return fn(
            new_loc,
            cds_frames=new_frames,
            sequence_guid=self.sequence_guid,
            sequence_name=self.sequence_name,
            protein_id=self.protein_id,
            product=self.product,
            qualifiers=self._export_qualifiers_to_list(),
            guid=None,  # generate a new Interval GUID based on updated data
        )
