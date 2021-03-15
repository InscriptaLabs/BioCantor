from enum import Enum


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
