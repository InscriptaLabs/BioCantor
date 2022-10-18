import pytest

from inscripta.biocantor.gene.codon import Codon, TranslationTable


class TestCodon:
    def test__init__(self):
        # case-insensitive
        obs = Codon("tCg")
        assert obs.value == "TCG"
        assert obs.name == "TCG"

    @pytest.mark.parametrize(
        "codon_str, exp_err",
        [
            ("A", ValueError),
            ("AG", ValueError),
            ("AGTR", ValueError),
            ("TCCTA", ValueError),
            ("GC-", ValueError),
            ("FRT", ValueError),
        ],
    )
    def test__init__errors(self, codon_str, exp_err):
        with pytest.raises(exp_err):
            Codon(codon_str)

    def test_singletons(self):
        obs1 = Codon("AAA")
        obs2 = Codon("AAA")
        obs3 = Codon("AAT")
        assert obs1 is obs2
        assert obs2 is not obs3

    def test__repr__(self):
        assert repr(Codon("GcT")) == "<Codon.GCT: GCT>"

    def test__str__(self):
        assert str(Codon("Agt")) == "AGT"

    def test__hash__(self):
        assert hash(Codon("GgC")) == hash("GGC")

    @pytest.mark.parametrize(
        "codon,expected",
        [
            (Codon("AAA"), "K"),
            (Codon("TGA"), "*"),
            # Wobble base N but translatable
            (Codon("GCN"), "A"),
            # Untranslatable
            (Codon("NGC"), "X"),
            (Codon("NNN"), "X"),
        ],
    )
    def test_translate(self, codon, expected):
        assert codon.translate() is expected

    @pytest.mark.parametrize(
        "codon, include_self, expected",
        [
            (Codon("GGT"), True, [Codon("GGT"), Codon("GGC"), Codon("GGA"), Codon("GGG")]),
            (Codon("GGT"), False, [Codon("GGC"), Codon("GGA"), Codon("GGG")]),
            # Don't keep the N if translatable codon, even if include_self is True
            (Codon("GCN"), True, [Codon("GCA"), Codon("GCG"), Codon("GCT"), Codon("GCC")]),
            (Codon("GCN"), False, [Codon("GCA"), Codon("GCG"), Codon("GCT"), Codon("GCC")]),
            # Keep self if untranslatable
            (Codon("NNN"), True, [Codon("NNN")]),
            (Codon("NNN"), False, []),
        ],
    )
    def test_synonymous_codons(self, codon, include_self, expected):
        assert set(codon.synonymous_codons(include_self)) == set(expected)

    @pytest.mark.parametrize(
        "codon_str, exp",
        [
            ("ATG", False),
            ("TgA", True),
            ("Tag", True),
            ("TAA", True),
            ("TAN", False),
        ],
    )
    def test_is_stop_codon(self, codon_str, exp):
        assert Codon(codon_str).is_stop_codon is exp

    @pytest.mark.parametrize(
        "codon_str, exp",
        [
            ("ATG", True),
            ("TgA", True),
            ("TNg", False),
            ("WRG", False),
            ("ATS", False),
        ],
    )
    def test_is_strict_codon(self, codon_str, exp):
        assert Codon(codon_str).is_strict_codon is exp

    @pytest.mark.parametrize(
        "codon_str, exp",
        [
            ("ATG", True),
            ("TgA", False),
        ],
    )
    def test_is_cannonical_start_codon(self, codon_str, exp):
        assert Codon(codon_str).is_canonical_start_codon is exp

    @pytest.mark.parametrize(
        "codon_str, table, exp",
        [
            ("ATG", TranslationTable.DEFAULT, True),
            ("ATG", TranslationTable.STANDARD, True),
            ("ATG", TranslationTable.PROKARYOTE, True),
            ("TTG", TranslationTable.DEFAULT, False),
            ("TTG", TranslationTable.STANDARD, True),
            ("TTG", TranslationTable.PROKARYOTE, True),
            ("ATA", TranslationTable.DEFAULT, False),
            ("ATA", TranslationTable.STANDARD, False),
            ("ATA", TranslationTable.PROKARYOTE, True),
            ("GGG", TranslationTable.DEFAULT, False),
            ("GGG", TranslationTable.STANDARD, False),
            ("GGG", TranslationTable.PROKARYOTE, False),
        ],
    )
    def test_is_start_codon_in_specific_translation_table(self, codon_str, table, exp):
        assert Codon(codon_str).is_start_codon_in_specific_translation_table(table) is exp
