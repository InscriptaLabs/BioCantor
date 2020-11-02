from enum import Enum


class Alphabet(Enum):
    NT_STRICT = "ACGT"
    NT_EXTENDED = "ATUCGNWSMKRYBDHV"
    NT_STRICT_GAPPED = "ACGT-"
    NT_EXTENDED_GAPPED = "ATUCGNWSMKRYBDHV-"
    AA = "GALMFWKQESPVICYHRNDT*"
    GENERIC = "ABCDEFGHIJKLMNOPQRSTUVWXYZ-"

    def is_nucleotide_alphabet(self):
        if self in [
            Alphabet.NT_STRICT,
            Alphabet.NT_EXTENDED,
            Alphabet.NT_STRICT_GAPPED,
            Alphabet.NT_EXTENDED_GAPPED,
        ]:
            return True
        if self in [Alphabet.AA, Alphabet.GENERIC]:
            return False
        raise NotImplementedError("Not implemented for alphabet {}".format(self))


ALPHABET_TO_NUCLEOTIDE_COMPLEMENT = {
    Alphabet.NT_STRICT: {"A": "T", "a": "t", "C": "G", "c": "g", "G": "C", "g": "c", "T": "A", "t": "a"},
    Alphabet.NT_EXTENDED: {
        "A": "T",
        "a": "t",
        "T": "A",
        "t": "a",
        "U": "A",
        "u": "a",
        "G": "C",
        "g": "c",
        "C": "G",
        "c": "g",
        "Y": "R",
        "y": "r",
        "R": "Y",
        "r": "y",
        "S": "S",
        "s": "s",
        "W": "W",
        "w": "w",
        "K": "M",
        "k": "m",
        "M": "K",
        "m": "k",
        "B": "V",
        "b": "v",
        "D": "H",
        "d": "h",
        "H": "D",
        "h": "d",
        "V": "B",
        "v": "b",
        "N": "N",
        "n": "n",
    },
    Alphabet.NT_STRICT_GAPPED: {
        "A": "T",
        "a": "t",
        "C": "G",
        "c": "g",
        "G": "C",
        "g": "c",
        "T": "A",
        "t": "a",
        "-": "-",
    },
    Alphabet.NT_EXTENDED_GAPPED: {
        "A": "T",
        "a": "t",
        "T": "A",
        "t": "a",
        "U": "A",
        "u": "a",
        "G": "C",
        "g": "c",
        "C": "G",
        "c": "g",
        "Y": "R",
        "y": "r",
        "R": "Y",
        "r": "y",
        "S": "S",
        "s": "s",
        "W": "W",
        "w": "w",
        "K": "M",
        "k": "m",
        "M": "K",
        "m": "k",
        "B": "V",
        "b": "v",
        "D": "H",
        "d": "h",
        "H": "D",
        "h": "d",
        "V": "B",
        "v": "b",
        "N": "N",
        "n": "n",
        "-": "-",
    },
}
