from enum import IntEnum
from typing import List, Optional, Union

from inscripta.biocantor.constants import gencode, extended_gencode, aacodons


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


class Codon:
    """Enum-like class for dealing with Codons"""

    __slots__ = ["_val"]
    _singletons_ = {}

    def __new__(cls, codon: Union[str, "Sequence"]):  # noqa: F821
        clean_codon = str(codon).upper()
        if clean_codon in cls._singletons_:
            return cls._singletons_[clean_codon]
        instance = super().__new__(cls)
        cls._singletons_[clean_codon] = instance
        return instance

    def __init__(self, codon: Union[str, "Sequence"]):  # noqa: F821
        self._val = str(codon).upper()
        if len(self._val) != 3:
            raise ValueError("Codon not a multiple of 3")
        if self._val.strip("ATUCGNWSMKRYBDHV") != "":
            raise ValueError(f"Unknown, non-nucleotide bases given: '{self._val}'")

    def __eq__(self, other):
        return other is self

    def __repr__(self):
        return f"<Codon.{self._val}: {self._val}>"

    def __str__(self):
        return self._val

    def __hash__(self):
        return hash(self._val)

    @property
    def value(self):
        return self._val

    @property
    def name(self):
        return self._val

    def translate(self) -> str:
        """Returns string symbol of translated amino acid"""
        try:
            return gencode[self._val]
        except KeyError:
            return extended_gencode[self._val] if self._val in extended_gencode else "X"

    def synonymous_codons(self, include_self=False) -> List["Codon"]:
        """Returns list of synonymous codons

        Parameters
        ----------
        include_self
            Include this codon in returned list
        """
        aa = self.translate()
        if aa == "X":
            return [self] if include_self else []

        return [Codon(codon_str) for codon_str in aacodons[aa] if include_self or codon_str != self._val]

    @property
    def is_stop_codon(self) -> bool:
        return self._val in aacodons["*"]

    @property
    def is_strict_codon(self) -> bool:
        return self._val in gencode

    @property
    def is_canonical_start_codon(self) -> bool:
        """Is this a canonical start codon?"""
        return self._val == "ATG"

    def is_start_codon_in_specific_translation_table(
        self, translation_table: Optional[TranslationTable] = TranslationTable.DEFAULT
    ) -> bool:
        """Prokaryotes have a wider pool of valid start codons than eukaryotes who use the canonical ATG only"""
        return self in START_CODONS_BY_TRANSLATION_TABLE[translation_table]


# NCBI maintains a set of translation tables with integer indices
# Index 1 is standard (vertebrate), while index 11 is prokaryotes
START_CODONS_BY_TRANSLATION_TABLE = {
    TranslationTable.DEFAULT: frozenset({Codon("ATG")}),
    TranslationTable.STANDARD: frozenset({Codon("ATG"), Codon("TTG"), Codon("CTG")}),
    TranslationTable.PROKARYOTE: frozenset(
        {Codon("ATG"), Codon("TTG"), Codon("CTG"), Codon("ATT"), Codon("ATC"), Codon("ATA"), Codon("GTG")}
    ),
}
