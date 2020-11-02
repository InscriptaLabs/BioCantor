from enum import Enum
from typing import List

from inscripta.biocantor.constants import gencode, aacodons


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
