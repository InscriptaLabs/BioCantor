import pytest

from inscripta.biocantor.gene.codon import Codon


class TestCodon:
    @pytest.mark.parametrize(
        "codon,expected",
        [
            (Codon.AAA, "K"),
            (Codon.TGA, "*"),
        ],
    )
    def test_translate(self, codon, expected):
        assert codon.translate() is expected

    @pytest.mark.parametrize(
        "codon, include_self, expected",
        [
            (Codon.GGT, True, [Codon.GGT, Codon.GGC, Codon.GGA, Codon.GGG]),
            (Codon.GGT, False, [Codon.GGC, Codon.GGA, Codon.GGG]),
        ],
    )
    def test_synonymous_codons(self, codon, include_self, expected):
        assert set(codon.synonymous_codons(include_self)) == set(expected)
