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

    def __eq__(self, other) -> bool:
        """Check equality through singleton comparison"""
        return other is self

    def __repr__(self) -> str:
        """Representation string of codon, equal to what would come out of an Enum object"""
        return f"<Codon.{self._val}: {self._val}>"

    def __str__(self) -> str:
        """The codon as a string"""
        return self._val

    def __hash__(self) -> int:
        """Hash of the codon string value"""
        return hash(self._val)

    @property
    def value(self) -> str:
        """The string value of the codon, to maintain enum-like functionality"""
        return self._val

    @property
    def name(self) -> str:
        """The string value of the codon, to maintain enum-like functionality"""
        return self._val

    def translate(self, strict: bool = True) -> str:
        """Returns string symbol of translated amino acid

        Parameters
        ----------
        strict
            Whether to only use strict ATGC codon translations, or allow translation using extended IUPAC sequence.
            Untranslatable amino acids, whether through strict or simply unknown, will be represented with an "X"
            Default True (strict ATGC only translation)
        """
        try:
            return gencode[self._val]
        except KeyError:
            if not strict and self._val in extended_gencode:
                return extended_gencode[self._val]
            else:
                return "X"

    def synonymous_codons(self, include_self: bool = False) -> List["Codon"]:
        """Returns list of synonymous codons

        Parameters
        ----------
        include_self
            Include this codon in returned list

        Returns
        -------
        List[Codon]
            The synonymous codons for the current codon.

        Notes
        -----
        For extended IUPAC codons, this function will still try and translate the codon in an attempt to get
        synonymous codons. If it can, the function will send out all strict ATGC codons for that amino acid.
        If it can't, it will send out either no codons or just self, if include_self is True.
        """
        aa = self.translate(strict=False)
        if aa == "X":
            return [self] if include_self else []

        return [Codon(codon_str) for codon_str in aacodons[aa] if include_self or codon_str != self._val]

    @property
    def is_stop_codon(self) -> bool:
        """Whether the codon is a stop codon or not

        Returns
        -------
        bool
            Whether the codon is a stop codon (True) or not (False)
        """
        return self._val in aacodons["*"]

    @property
    def is_strict_codon(self) -> bool:
        """Whether the codon is a strict ATGC containing codon or not

        Returns
        -------
        bool
            Whether the codon is a strict ATGC codon (True) or not (False)
        """
        return self._val in gencode

    @property
    def is_canonical_start_codon(self) -> bool:
        """Whether the codon is a canonical start codon or not

        Returns
        -------
        bool
            Whether the codon is a canonical start codon (True) or not (False)
        """
        return self._val == "ATG"

    def is_start_codon_in_specific_translation_table(
        self, translation_table: Optional[TranslationTable] = TranslationTable.DEFAULT
    ) -> bool:
        """Whether the codon is a start codon according to the given translation_table

        Parameters
        ----------
        translation_table
            What translation table to check for start codon validation.

        Returns
        -------
        bool
            Whether the codon is a start codon (True) or not (False)

        Notes
        -----
        Prokaryotes have a wider pool of valid start codons than eukaryotes who use the canonical ATG only
        """
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
