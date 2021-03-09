from enum import Enum, IntEnum
from typing import List, Optional

from inscripta.biocantor.constants import gencode, aacodons


class TranslationTable(IntEnum):
    """
    NCBI maintains a set of translation tables:

    https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG11

    For most purposes, the only two tables that are relevant are:

    * Table 1: The Standard Code (vertebrates)
    * Table 11: The Bacterial, Archaeal and Plant Plastid Code

    The only difference between these two tables are the choice of initiator codons. While the standard start codon
    is `ATG`, in certain eukaryotes it is possible for `TTG` and `CTG` to be start codons.

    In table 11, it is possible for any of the following codons to be considered initiator codons:

    .. code-block::
        ATG, TTG, CTG, ATT, ATC, ATA, GTG

    These alternative initiator codons do see usage regularly, with approximately 3% of E. coli using TTG.

    For this enumeration, the above two tables are included, as well as the default ATG-only option. This is because
    for many uses, ATG is the only one worth considering.
    """

    DEFAULT = 0
    STANDARD = 1
    PROKARYOTE = 11


class Codon(Enum):

    ATA = "ATA"
    ATC = "ATC"
    ATT = "ATT"
    ATG = "ATG"
    ACA = "ACA"
    ACC = "ACC"
    ACG = "ACG"
    ACT = "ACT"
    AAC = "AAC"
    AAT = "AAT"
    AAA = "AAA"
    AAG = "AAG"
    AGC = "AGC"
    AGT = "AGT"
    AGA = "AGA"
    AGG = "AGG"
    CTA = "CTA"
    CTC = "CTC"
    CTG = "CTG"
    CTT = "CTT"
    CCA = "CCA"
    CCC = "CCC"
    CCG = "CCG"
    CCT = "CCT"
    CAC = "CAC"
    CAT = "CAT"
    CAA = "CAA"
    CAG = "CAG"
    CGA = "CGA"
    CGC = "CGC"
    CGG = "CGG"
    CGT = "CGT"
    GTA = "GTA"
    GTC = "GTC"
    GTG = "GTG"
    GTT = "GTT"
    GCA = "GCA"
    GCC = "GCC"
    GCG = "GCG"
    GCT = "GCT"
    GAC = "GAC"
    GAT = "GAT"
    GAA = "GAA"
    GAG = "GAG"
    GGA = "GGA"
    GGC = "GGC"
    GGG = "GGG"
    GGT = "GGT"
    TCA = "TCA"
    TCC = "TCC"
    TCG = "TCG"
    TCT = "TCT"
    TTC = "TTC"
    TTT = "TTT"
    TTA = "TTA"
    TTG = "TTG"
    TAC = "TAC"
    TAT = "TAT"
    TAA = "TAA"
    TAG = "TAG"
    TGC = "TGC"
    TGT = "TGT"
    TGA = "TGA"
    TGG = "TGG"

    def translate(self) -> str:
        """Returns string symbol of translated amino acid"""
        return gencode[self.value]

    def synonymous_codons(self, include_self=False) -> List["Codon"]:
        """Returns list of synonymous codons

        Parameters
        ----------
        include_self
            Include this codon in returned list
        """
        codons = [Codon(codon_str) for codon_str in aacodons[self.translate()]]
        if not include_self:
            codons.remove(self)
        return codons

    @property
    def is_stop_codon(self) -> bool:
        return self.value in aacodons["*"]

    @property
    def is_canonical_start_codon(self) -> bool:
        """Is this a canonical start codon?"""
        return self == Codon.ATG

    def is_start_codon_in_specific_translation_table(
        self, translation_table: Optional[TranslationTable] = TranslationTable.DEFAULT
    ) -> bool:
        """Prokaryotes have a wider pool of valid start codons than eukaryotes who use the canonical ATG only"""
        return self in START_CODONS_BY_TRANSLATION_TABLE[translation_table]


# NCBI maintains a set of translation tables with integer indices
# Index 1 is standard (vertebrate), while index 11 is prokaryotes
START_CODONS_BY_TRANSLATION_TABLE = {
    TranslationTable.DEFAULT: {Codon.ATG},
    TranslationTable.STANDARD: {Codon.ATG, Codon.TTG, Codon.CTG},
    TranslationTable.PROKARYOTE: {Codon.ATG, Codon.TTG, Codon.CTG, Codon.ATT, Codon.ATC, Codon.ATA, Codon.GTG},
}
