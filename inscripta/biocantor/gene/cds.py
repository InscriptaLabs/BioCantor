from enum import Enum
from methodtools import lru_cache
from typing import Iterable, List, Union, Optional

from inscripta.biocantor.exc import LocationException
from inscripta.biocantor.gene.codon import Codon, TranslationTable
from inscripta.biocantor.location.location import Location, Strand
from inscripta.biocantor.location.location_impl import SingleInterval, CompoundInterval
from inscripta.biocantor.sequence import Sequence
from inscripta.biocantor.sequence.alphabet import Alphabet


class CDSPhase(Enum):
    """
    It is important not to confuse Phase with Frame. From the GFF3 specification:

    The phase is one of the integers 0, 1, or 2, indicating the number of bases forward from the start of the
    current CDS feature the next codon begins. A phase of "0" indicates that a codon begins on the first nucleotide
    of the CDS feature (i.e. 0 bases forward), a phase of "1" indicates that the codon begins at the second nucleotide
    of this CDS feature and a phase of "2" indicates that the codon begins at the third nucleotide of this region.

    """

    NONE = -1
    ZERO = 0
    ONE = 1
    TWO = 2

    @staticmethod
    def from_int(value: int) -> "CDSPhase":
        return CDSPhase(value)  # Raises ValueError for invalid int

    def to_frame(self) -> "CDSFrame":
        """
        https://github.com/ucscGenomeBrowser/kent/blob/022eb4f62a0af16526ca1bebcd9e68bd456265dc/src/inc/gff3.h#L281-L293
        """
        mapping = {0: 0, 2: 1, 1: 2, -1: -1}
        return CDSFrame(mapping[self.value])

    def to_gff(self) -> str:
        """In GFF format, Phase is represented with a period for NONE"""
        if self == CDSPhase.NONE:
            return "."
        return str(self.value)


class CDSFrame(Enum):
    """
    From the GFF3 specification:

    Frame is generally calculated as a value for a given base relative to the start of the complete
    open reading frame (ORF) or the codon (e.g. modulo 3) while CDS phase describes the start of the next codon
    relative to a given CDS feature.

    Frame is easier to work with computationally because for any given position it can be calculated by the distance
    from transcription start % 3. Care must be taken for strand -- if the CDS is on the minus strand, then this
    calculation instead becomes distance from transcription stop % 3.
    """

    NONE = -1
    ZERO = 0
    ONE = 1
    TWO = 2

    @staticmethod
    def from_int(value: int) -> "CDSFrame":
        return CDSFrame(value)

    def shift(self, shift: int) -> "CDSFrame":
        if self is CDSFrame.NONE:
            return self
        if shift > 0:
            return CDSFrame.from_int((self.value + shift) % 3)
        else:
            return CDSFrame.from_int((self.value - (shift - ((-shift) % 3))) % 3)

    def to_phase(self) -> "CDSPhase":
        """Converts frame to phase"""
        mapping = {0: 0, 1: 2, 2: 1, -1: -1}
        return CDSPhase(mapping[self.value])


class CDSInterval:
    """
    An wrapper for a Location that gives it a concept of Frames.
    Frame must be recorded individually on *every interval*.

    This is necessary to be able to encode programmed frame shifts and indel variation that
    would break frame, but that are likely just errors in assembly or annotation.
    """

    def __init__(self, location: Location, frames: List[Union[CDSFrame, CDSPhase]]):
        self.location = location
        if len(frames) != location.num_blocks:
            raise LocationException("Number of frame entries must match number of exons")
        # internally we work with Frame, but support Phase
        # this will make parsing GFF3 easier
        for i, frame in enumerate(frames):
            if isinstance(frame, CDSPhase):
                frame = frame.to_frame()
                frames[i] = frame
        self.frames = frames

    def __str__(self):
        frame_str = ", ".join([str(p) for p in self.frames])
        return f"CDS(({self.location}), ({frame_str})"

    def __repr__(self):
        return "<{}>".format(str(self))

    def __eq__(self, other):
        if type(other) is not CDSInterval:
            return False
        elif self.frames != other.frames:
            return False
        return self.location == other.location

    def __hash__(self):
        return hash((self.location, self.frames[0]))

    def __len__(self) -> int:
        return len(self.location)

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

    @property
    def strand(self) -> Strand:
        """Pass up the Strand of this CDS's Location"""
        return self.location.strand

    @property
    def start(self) -> int:
        """Pass up the start of this CDS's Location"""
        return self.location.start

    @property
    def end(self) -> int:
        """Pass up the end of this CDS's Location"""
        return self.location.end

    def frame_iter(self) -> Iterable[CDSFrame]:
        """Iterate over frames taking strand into account"""
        if self.location.strand == Strand.PLUS or self.location.strand == Strand.UNSTRANDED:
            yield from self.frames
        else:
            yield from reversed(self.frames)

    def exon_iter(self) -> Iterable[SingleInterval]:
        """Iterate over exons in transcription direction"""
        if self.location.strand == Strand.PLUS or self.location.strand == Strand.UNSTRANDED:
            yield from self.location.blocks
        else:
            yield from reversed(self.location.blocks)

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

        Any leading or trailing bases that are annotated as CDS but cannot form a full codon
        are excluded. Additionally, any internal codons that are incomplete are excluded.

        Incomplete internal codons are determined by comparing the CDSFrame of each exon
        as annotated, to the expected value of the CDSFrame. This allows for an annotation
        to model things like programmed frameshifts and indels that may be assembly errors.
        """
        return len(list(self.scan_codon_locations()))

    def scan_codons(self, truncate_at_in_frame_stop: Optional[bool] = False) -> Iterable[Codon]:
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

    def scan_codon_locations(self) -> Iterable[Location]:
        """
        Returns an iterator over codon locations.

        Any leading or trailing bases that are annotated as CDS but cannot form a full codon
        are excluded. Additionally, any internal codons that are incomplete are excluded.

        Incomplete internal codons are determined by comparing the CDSFrame of each exon
        as annotated, to the expected value of the CDSFrame. This allows for an annotation
        to model things like programmed frameshifts and indels that may be assembly errors.
        """
        next_frame = CDSFrame.ZERO
        cleaned_rel_starts = []
        cleaned_rel_ends = []
        for exon, frame in zip(self.exon_iter(), self.frame_iter()):
            start_to_rel = self.location.parent_to_relative_pos(exon.start)
            end_to_rel_inclusive = self.location.parent_to_relative_pos(exon.end - 1)
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
            self.location.relative_interval_to_parent_location(cleaned_rel_starts[i], cleaned_rel_ends[i], Strand.PLUS)
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
        aa_seq_str = "".join([codon.translate() for codon in self.scan_codons(truncate_at_in_frame_stop)])
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

        ```
        CompoundInterval([0, 7, 12], [5, 11, 18], Strand.PLUS)
        ```

        ```
        Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
        Sequence:   A A A C A A A A G G G  T  A  C  C  C  A  A  A  A  A  A
        Exons:      A A A C A     A G G G     A  C  C  C  A  A
        Zero Frame: 0 1 2 0 1     2 0 1 2     0  1  2  0  1  2
        One Frame:  - 0 1 2 0     1 2 0 1     2  0  1  2  0  1
        Two Frame:  - - 0 1 2     0 1 2 0     1  2  0  1  2  0
        ```

        In the non-zero case, the ``[0, 1, 2]`` cycle is offset by 1 or 2 bases.

        So, for this test case we expect the frames to be:

        ```
        Zero Frame: [0, 2, 0]
        One Frame:  [1, 1, 2]
        Two Frame:  [2, 0, 1]
        ```

        2. Minus strand:

        ```
        CompoundInterval([0, 7, 12], [5, 11, 18], Strand.MINUS)
        ```

        ```
        Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
        Sequence:   A A A C A A A A G G G  T  A  C  C  C  A  A  A  A  A  A
        Exons:      A A A C A     A G G G     A  C  C  C  A  A
        Zero Frame: 2 1 0 2 1     0 2 1 0     2  1  0  2  1  0
        One Frame:  1 0 2 1 0     2 1 0 2     1  0  2  1  0  -
        Two Frame:  0 2 1 0 2     1 0 2 1     0  2  1  0  -  -
        ```

        Now, for negative strand CDS intervals, the frame list is still in plus strand orientation.

        So, for this test case we expect the frames to be:

        ```
        Zero Frame: [1, 0, 0]
        One Frame:  [0, 2, 1]
        Two Frame:  [2, 1, 2]
        ```


        Args:
            location: A interval of the CDS.
            starting_frame: Frame to start iteration with. If ``codon_start`` was the source of this value,
                then you would subtract one before converting to :class:`CDSFrame`.

        Returns:
            A list of :class:`CDSFrame` that could be combined with the input Location to build a :class:`CDSInterval`.
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
        new_loc = self.location.optimize_blocks()
        first_frame = next(self.frame_iter())
        frames = CDSInterval.construct_frames_from_location(new_loc, first_frame)
        return CDSInterval(new_loc, frames)

    def optimize_and_combine_blocks(self) -> "CDSInterval":
        """
        Combine the blocks of this CDS interval, including removing overlapping blocks.

        Once this operation is performed, internal frameshifts modeled by 0bp gaps will be lost, as well as programmed
        frameshifts modeled by overlapping blocks. The resulting translations will be out of frame downstream.

        Returns:
            A new :class:`CDSInterval` that has been merged.
        """
        new_loc = self.location.optimize_and_combine_blocks()
        first_frame = next(self.frame_iter())
        frames = CDSInterval.construct_frames_from_location(new_loc, first_frame)
        return CDSInterval(new_loc, frames)
