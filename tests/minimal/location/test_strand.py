import pytest

from biocantor.exc import InvalidStrandException
from biocantor.location.strand import Strand


class TestStrand:
    @pytest.mark.parametrize(
        "strand,string",
        [
            (Strand.PLUS, "+"),
            (Strand.MINUS, "-"),
            (Strand.UNSTRANDED, "."),
        ],
    )
    def test_str(self, strand, string):
        assert str(strand) == string

    def test_from_symbol_valid(self):
        assert Strand.from_symbol("+") == Strand.PLUS
        assert Strand.from_symbol("-") == Strand.MINUS

    def test_from_symbol_invalid(self):
        with pytest.raises(ValueError):
            Strand.from_symbol("x")

    def test_to_symbol(self):
        assert Strand.PLUS.to_symbol() == "+"
        assert Strand.MINUS.to_symbol() == "-"
        assert Strand.UNSTRANDED.to_symbol() == "."

    def test_from_int_valid(self):
        assert Strand.from_int(1) == Strand.PLUS
        assert Strand.from_int(-1) == Strand.MINUS
        assert Strand.from_int(0) == Strand.UNSTRANDED

    def test_from_int_invalid(self):
        with pytest.raises(ValueError):
            Strand.from_int(2)

    def test_equals(self):
        s = {Strand.PLUS, Strand.from_symbol("+")}
        assert len(s) == 1

    def test_order(self):
        assert sorted([Strand.MINUS, Strand.UNSTRANDED, Strand.PLUS]) == [
            Strand.PLUS,
            Strand.MINUS,
            Strand.UNSTRANDED,
        ]

    def test_reverse(self):
        assert Strand.PLUS.reverse() == Strand.MINUS
        assert Strand.MINUS.reverse() == Strand.PLUS
        assert Strand.UNSTRANDED.reverse() == Strand.UNSTRANDED

    def test_relative_to(self):
        assert Strand.PLUS.relative_to(Strand.PLUS) == Strand.PLUS
        assert Strand.PLUS.relative_to(Strand.MINUS) == Strand.MINUS
        assert Strand.PLUS.relative_to(Strand.UNSTRANDED) == Strand.UNSTRANDED
        assert Strand.MINUS.relative_to(Strand.PLUS) == Strand.MINUS
        assert Strand.MINUS.relative_to(Strand.MINUS) == Strand.PLUS
        assert Strand.MINUS.relative_to(Strand.UNSTRANDED) == Strand.UNSTRANDED
        assert Strand.UNSTRANDED.relative_to(Strand.PLUS) == Strand.UNSTRANDED
        assert Strand.UNSTRANDED.relative_to(Strand.MINUS) == Strand.UNSTRANDED
        assert Strand.UNSTRANDED.relative_to(Strand.UNSTRANDED) == Strand.UNSTRANDED

    def test_assert_directional(self):
        Strand.PLUS.assert_directional()
        Strand.MINUS.assert_directional()
        with pytest.raises(InvalidStrandException):
            Strand.UNSTRANDED.assert_directional()
