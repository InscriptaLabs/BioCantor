from itertools import count, zip_longest
from typing import Iterator, List, Union, Optional, Dict, Hashable, Any, Iterable, Set
from uuid import UUID

from methodtools import lru_cache

from inscripta.biocantor.exc import (
    InvalidCDSIntervalError,
    NoSuchAncestorException,
    LocationOverlapException,
    MismatchedFrameException,
)
from inscripta.biocantor.gene.cds_frame import CDSPhase, CDSFrame
from inscripta.biocantor.gene.codon import Codon, TranslationTable
from inscripta.biocantor.gene.interval import AbstractFeatureInterval, QualifierValue
from inscripta.biocantor.io.bed import RGB, BED12
from inscripta.biocantor.io.gff3.constants import GFF_SOURCE, NULL_COLUMN, BioCantorFeatureTypes, BioCantorQualifiers
from inscripta.biocantor.io.gff3.rows import GFFAttributes, GFFRow
from inscripta.biocantor.location.location import Location, Strand
from inscripta.biocantor.location.location_impl import SingleInterval, CompoundInterval
from inscripta.biocantor.parent.parent import Parent, SequenceType
from inscripta.biocantor.sequence import Sequence
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.util.hashing import digest_object


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
                self.frames,
                self.product,
                self.protein_id,
                self.qualifiers,
            )
        else:
            self.guid = guid

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

    @lru_cache(maxsize=1)
    @property
    def chunk_relative_frames(self) -> List[CDSFrame]:
        """
        It may be the case that the chunk relative location of this CDSInterval object is a subset
        of the full chromosomal location. In this case, the frames list needs to be appropriately
        subsetted to the correct set of frame entries.

        Returns:
            A list of CDSFrame entries that overlap the chunk relative location.
        """
        frames = []
        for genomic_start, genomic_end, frame in zip(self._genomic_starts, self._genomic_ends, self.frames):
            genomic_exon = SingleInterval(
                genomic_start, genomic_end, Strand.PLUS, parent=self.chromosome_location.parent
            )
            # chromosome location has overlapping blocks merged, so that the intersection always has one block
            # this is OK to do here since the original genomic intervals retain the overlapping information
            if isinstance(self.chromosome_location, SingleInterval):
                chrom_loc = self.chromosome_location
            elif isinstance(self.chromosome_location, CompoundInterval):
                chrom_loc = self.chromosome_location.optimize_and_combine_blocks()
            else:
                return frames
            intersection = genomic_exon.intersection(chrom_loc, match_strand=False)

            if intersection.is_empty:
                continue
            elif intersection.num_blocks != 1:
                raise LocationOverlapException("Found overlapping blocks after block optimization")
            elif len(intersection) == len(genomic_exon):
                frames.append(frame)
            else:
                # shift frame forwards by the amount the exon was sliced from the left
                shift = intersection.start - genomic_exon.start
                new_frame = frame.shift(-shift)
                frames.append(new_frame)
        return frames

    def to_dict(self, chromosome_relative_coordinates: bool = True) -> Dict[str, Any]:
        """
        Convert this CDS to a dictionary representation. Note that the dictionary representation can only
        use CDSFrame, not CDSPhase.
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
        :meth:`CDSInterval.to_dict()`. The frames must be CDSFrame, not CDSPhase.
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
        present on a sequence chunk, which could be a sequence produced

        This location should
        be built by something like this:

        .. code-block:: python

            from inscripta.biocantor.io.parser import seq_chunk_to_parent
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

    def to_gff(
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterable[GFFRow]:

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
            attributes = GFFAttributes(
                id=f"{cds_guid}-{i}",
                qualifiers=qualifiers,
                name=self.protein_id,
                parent=parent,
            )
            row = GFFRow(
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

    def frame_iter(self, chunk_relative_frames: bool = True) -> Iterator[CDSFrame]:
        """Iterate over frames taking strand into account.

        If ``chunk_relative_frames`` is ``True``, then this iterator will only iterate over frames that
        overlap the relative chunk. These frames will potentially be reduced in quantity, and also shifted to handle
        exons that are now partial exons.
        """
        frames = self.frames if chunk_relative_frames is False else self.chunk_relative_frames
        if (
            self.chunk_relative_location.strand == Strand.PLUS
            or self.chunk_relative_location.strand == Strand.UNSTRANDED
        ):
            yield from frames
        else:
            yield from reversed(frames)

    def exon_iter(self) -> Iterator[SingleInterval]:
        """Iterate over exons in transcription direction"""
        if (
            self.chunk_relative_location.strand == Strand.PLUS
            or self.chunk_relative_location.strand == Strand.UNSTRANDED
        ):
            yield from self._location.blocks
        else:
            yield from reversed(self.chunk_relative_location.blocks)

    @lru_cache(maxsize=1)
    def extract_sequence(self) -> Sequence:
        """
        Returns a continuous CDS sequence that is in frame and always a multiple of 3.

        Any leading or trailing bases that are annotated as CDS but cannot form a full codon
        are removed. Additionally, any internal codons that are incomplete are removed.

        Incomplete internal codons are determined by comparing the CDSFrame of each exon
        as annotated, to the expected value of the CDSFrame. This allows for an annotation
        to model things like programmed frameshifts and indels that may be assembly errors.
        """
        codons = [str(codon_location.extract_sequence()) for codon_location in self.scan_codon_locations()]
        seq = "".join(codons)
        assert len(seq) % 3 == 0
        return Sequence(seq, Alphabet.NT_EXTENDED)

    @property
    def num_codons(self) -> int:
        """
        Returns the number of codons.

        NOTE: If this CDS is a subset of the original sequence, this number will represent the subset,
        not the original size!

        TODO: Make this more efficient.

        Any leading or trailing bases that are annotated as CDS but cannot form a full codon
        are excluded. Additionally, any internal codons that are incomplete are excluded.

        Incomplete internal codons are determined by comparing the CDSFrame of each exon
        as annotated, to the expected value of the CDSFrame. This allows for an annotation
        to model things like programmed frameshifts and indels that may be assembly errors.
        """
        return len(list(self.scan_codon_locations()))

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

    def scan_codon_locations(self) -> Iterator[Location]:
        """
        Returns an iterator over codon locations in *chunk relative* coordinates.

        TODO: Allow chromosome relative codon scanning. This is an issue because lifting to chromosome coordinates
            removes sequence information.

        Any leading or trailing bases that are annotated as CDS but cannot form a full codon
        are excluded. Additionally, any internal codons that are incomplete are excluded.

        Incomplete internal codons are determined by comparing the CDSFrame of each exon
        as annotated, to the expected value of the CDSFrame. This allows for an annotation
        to model things like programmed frameshifts and indels that may be assembly errors.
        """
        next_frame = CDSFrame.ZERO
        cleaned_rel_starts = []
        cleaned_rel_ends = []
        # zip_longest is used here to ensure that the two iterators are always actually in sync
        for exon, frame in zip_longest(self.exon_iter(), self.frame_iter()):

            if not exon or not frame:
                raise MismatchedFrameException("Frame iterator is not in sync with exon iterator")

            start_to_rel = self._location.parent_to_relative_pos(exon.start)
            end_to_rel_inclusive = self._location.parent_to_relative_pos(exon.end - 1)
            rel_start = min(start_to_rel, end_to_rel_inclusive)
            rel_end = max(start_to_rel, end_to_rel_inclusive) + 1
            if next_frame != frame:
                rel_start += frame.value
                # remove trailing codon from previous block
                shift = self._total_block_len(cleaned_rel_starts, cleaned_rel_ends) % 3
                if shift > 0:
                    cleaned_rel_ends[-1] = cleaned_rel_ends[-1] - shift
                # we are now inherently in frame
                next_frame = CDSFrame(CDSFrame.ZERO)
            # it may be the case that the removal of the trailing codon from the previous block entirely
            # eliminates that block. In that case, skip it.
            if rel_start >= rel_end:
                continue
            cleaned_rel_starts.append(rel_start)
            cleaned_rel_ends.append(rel_end)
            # this is what we expect the next frame to be, if no frameshift occurred
            next_frame = next_frame.shift(rel_end - rel_start)
        cleaned_blocks = [
            self._location.relative_interval_to_parent_location(cleaned_rel_starts[i], cleaned_rel_ends[i], Strand.PLUS)
            for i in range(len(cleaned_rel_starts))
        ]
        cleaned_location = CompoundInterval.from_single_intervals(cleaned_blocks)
        if len(cleaned_location) < 3:
            return
        yield from cleaned_location.scan_windows(3, 3, 0)

    @staticmethod
    def _total_block_len(starts: List[int], ends: List[int]) -> int:
        return sum([coords[1] - coords[0] for coords in zip(starts, ends)])

    @lru_cache(maxsize=2)
    def translate(self, truncate_at_in_frame_stop: Optional[bool] = False) -> Sequence:
        """
        Returns amino acid sequence of this CDS. If truncate_at_in_frame_stop is ``True``,
        this will stop at the first in-frame stop.
        """
        aa_seq_str = "".join((codon.translate() for codon in self.scan_codons(truncate_at_in_frame_stop)))
        return Sequence(aa_seq_str, Alphabet.AA)

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
        first_frame = next(self.frame_iter())
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
        new_loc = self.chunk_relative_location.optimize_and_combine_blocks()
        first_frame = next(self.frame_iter())
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
