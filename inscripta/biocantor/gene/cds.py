from enum import Enum
from typing import Iterator, List, Union, Optional, AnyStr

from inscripta.biocantor.exc import ValidationError, EmptyLocationException
from inscripta.biocantor.gene.codon import Codon
from inscripta.biocantor.location.location import Location, Strand
from inscripta.biocantor.location.location_impl import SingleInterval
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.sequence import Sequence


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

    def to_gff(self) -> AnyStr:
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
            raise ValidationError("Number of frame entries must match number of exons")
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

    def __len__(self) -> int:
        return len(self.location)

    @property
    def has_valid_stop(self) -> bool:
        """Does this CDS have a valid stop? Requires a sequence be associated."""
        seq = self.extract_sequence()
        c = Codon(seq[-3:].sequence)
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

    def intersect(self, location: Location) -> "CDSInterval":
        """
        Returns a new CDS representing the intersection of this CDS's location with the other location.
        Strand of the other location is ignored; returned CDS is on the same strand as this CDS.
        """
        intersection = self.location.intersection(location, match_strand=False)

        if intersection.is_empty:
            raise EmptyLocationException("Can't intersect disjoint intervals")

        frames = []
        if self.location.strand == Strand.PLUS or self.location.strand == Strand.UNSTRANDED:
            frame_shift = self.location.start - intersection.start
        else:
            frame_shift = intersection.end - self.location.end

        # adjust original frames by the new frame shift value
        for block, frame in zip(self.location.blocks, self.frames):
            if block.has_overlap(intersection):
                frames.append(frame.shift(frame_shift))

        return CDSInterval(intersection, frames)

    def frame_iter(self) -> Iterator[CDSFrame]:
        """Iterate over frames taking strand into account"""
        if self.location.strand == Strand.PLUS or self.location.strand == Strand.UNSTRANDED:
            yield from self.frames
        else:
            yield from reversed(self.frames)

    def exon_iter(self) -> Iterator[SingleInterval]:
        """Iterate over exons in transcription direction"""
        if self.location.strand == Strand.PLUS or self.location.strand == Strand.UNSTRANDED:
            yield from self.location.blocks
        else:
            yield from reversed(self.location.blocks)

    def extract_sequence(self) -> Sequence:
        """
        Returns a continuous CDS sequence that is in frame and always a multiple of 3.

        Any leading or trailing bases that are annotated as CDS but cannot form a full codon
        are removed. Additionally, any internal codons that are incomplete are removed.

        Incomplete internal codons are determined by comparing the CDSFrame of each exon
        as annotated, to the expected value of the CDSFrame. This allows for an annotation
        to model things like programmed frameshifts and indels that may be assembly errors.
        """
        seq = []
        # keeps track of what we expect the next frame to be
        next_frame = CDSFrame(CDSFrame.ZERO)
        for exon, frame in zip(self.exon_iter(), self.frame_iter()):
            s = list(exon.extract_sequence())
            if next_frame != frame:
                s = s[frame.value :]
                # remove trailing codon from previous sequence
                shift = len(seq) % 3
                if shift > 0:
                    seq = seq[:-shift]
                # we are now inherently in frame
                next_frame = CDSFrame(CDSFrame.ZERO)
            seq.extend(s)
            # this is what we expect the next frame to be, if no frameshift occurred
            next_frame = next_frame.shift(len(s))
        # we may not have a complete 3' end
        extra = len(seq) % 3
        if extra > 0:
            seq = "".join(str(s) for s in seq[:-extra])
        else:
            seq = "".join(str(s) for s in seq)
        assert len(seq) % 3 == 0
        return Sequence(seq, Alphabet.NT_EXTENDED)

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

    def translate(self, truncate_at_in_frame_stop: Optional[bool] = False) -> Sequence:
        """
        Returns amino acid sequence of this CDS. If truncate_at_in_frame_stop is True,
        this will stop at the first in-frame stop
        """
        aa_seq_str = "".join([codon.translate() for codon in self.scan_codons(truncate_at_in_frame_stop)])
        return Sequence(aa_seq_str, Alphabet.AA)
