import pytest

from inscripta.biocantor.exc import EmptyLocationException
from inscripta.biocantor.gene.cds import CDSPhase, CDSFrame, CDSInterval
from inscripta.biocantor.gene.codon import Codon
from inscripta.biocantor.location.location_impl import CompoundInterval, SingleInterval
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.sequence import Sequence


class TestCDSPhase:
    def test_from_int(self):
        assert CDSPhase.from_int(2) == CDSPhase.TWO
        with pytest.raises(ValueError):
            CDSPhase.from_int(3)

    @pytest.mark.parametrize(
        "phase,expected", [(CDSPhase.ONE, "1"), (CDSPhase.TWO, "2"), (CDSPhase.ZERO, "0"), (CDSPhase.NONE, ".")]
    )
    def test_to_gff(self, phase, expected):
        assert phase.to_gff() == expected


class TestCDSFrame:
    def test_from_int(self):
        assert CDSFrame.from_int(2) == CDSFrame.TWO
        with pytest.raises(ValueError):
            CDSFrame.from_int(3)

    @pytest.mark.parametrize(
        "frame,shift,expected",
        [
            (CDSFrame.ZERO, 0, CDSFrame.ZERO),
            (CDSFrame.ZERO, -1, CDSFrame.TWO),
            (CDSFrame.ZERO, 3, CDSFrame.ZERO),
            (CDSFrame.ONE, 1, CDSFrame.TWO),
            (CDSFrame.TWO, 5, CDSFrame.ONE),
        ],
    )
    def test_shift(self, frame, shift, expected):
        assert frame.shift(shift) == expected


class TestCDS:
    @pytest.mark.parametrize(
        "cds,location,expected",
        [
            (  # remove the first 2bp and last 2bp
                CDSInterval(SingleInterval(0, 10, Strand.PLUS), [CDSFrame.ZERO]),
                SingleInterval(2, 8, Strand.MINUS),
                CDSInterval(SingleInterval(2, 8, Strand.PLUS), [CDSFrame.ONE]),
            ),
            (  # remove the first 2bp and last 2bp; CDS on minus strand
                CDSInterval(SingleInterval(0, 10, Strand.MINUS), [CDSFrame.ZERO]),
                SingleInterval(2, 8, Strand.MINUS),
                CDSInterval(SingleInterval(2, 8, Strand.MINUS), [CDSFrame.ONE]),
            ),
            (  # strand doesn't matter
                CDSInterval(SingleInterval(5, 10, Strand.MINUS), [CDSFrame.TWO]),
                SingleInterval(0, 10, Strand.UNSTRANDED),
                CDSInterval(SingleInterval(5, 10, Strand.MINUS), [CDSFrame.TWO]),
            ),
            (  # slice off the first base
                CDSInterval(
                    CompoundInterval([4, 7, 12], [6, 10, 15], Strand.PLUS),
                    [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.TWO],
                ),
                SingleInterval(5, 15, Strand.UNSTRANDED),
                CDSInterval(
                    CompoundInterval([5, 7, 12], [6, 10, 15], Strand.PLUS),
                    [CDSFrame.TWO, CDSFrame.ONE, CDSFrame.ONE],
                ),
            ),
            (  # slice off the first exon
                CDSInterval(
                    CompoundInterval([4, 7, 12], [6, 10, 15], Strand.PLUS),
                    [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.TWO],
                ),
                SingleInterval(7, 15, Strand.UNSTRANDED),
                CDSInterval(
                    CompoundInterval([7, 12], [10, 15], Strand.PLUS),
                    [CDSFrame.TWO, CDSFrame.TWO],
                ),
            ),
            (  # slice off the last exon + 1bp
                CDSInterval(
                    CompoundInterval([4, 7, 12], [6, 10, 15], Strand.PLUS),
                    [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.TWO],
                ),
                SingleInterval(0, 8, Strand.UNSTRANDED),
                CDSInterval(
                    CompoundInterval([4, 7], [6, 8], Strand.PLUS),
                    [CDSFrame.ZERO, CDSFrame.TWO],
                ),
            ),
        ],
    )
    def test_intersect(self, cds, location, expected):
        assert cds.intersect(location) == expected

    @pytest.mark.parametrize(
        "cds,location,expected_exception",
        [
            (
                CDSInterval(SingleInterval(5, 10, Strand.PLUS), [CDSFrame.TWO]),
                SingleInterval(10, 15, Strand.PLUS),
                EmptyLocationException,
            ),
        ],
    )
    def test_intersect_error(self, cds, location, expected_exception):
        with pytest.raises(expected_exception):
            cds.intersect(location)

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # Contiguous CDS, plus strand, frame=0
            (
                CDSInterval(
                    SingleInterval(0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", Alphabet.NT_STRICT)),
                    [CDSFrame.ZERO],
                ),
                [Codon.ATA, Codon.CGA, Codon.TCA],
            ),
            # Discontiguous CDS, plus strand, frame=1, codons don't reach end of CDS
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", Alphabet.NT_STRICT),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                [Codon.CAG, Codon.GGA, Codon.CCC],  # QGP
            ),
            # Discontiguous CDS, plus strand, frame=1, 1bp deletion at start of exon 2
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [5, 16],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGACCCAAAAAA", Alphabet.NT_STRICT),
                    ),
                    [CDSFrame.ONE, CDSFrame.ZERO],
                ),
                [Codon.GGA, Codon.CCC],  # GP
            ),
            # Discontiguous CDS, plus strand, frame=1, 1bp insertion inside exon 2 (programmed frameshift)
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8, 12],
                        [5, 11, 18],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", Alphabet.NT_STRICT),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO, CDSFrame.TWO],  # QGP
                ),
                [Codon.CAG, Codon.GGA, Codon.CCC],  # QGP
            ),
            # Discontiguous CDS, minus strand, frame=2
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_STRICT),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                [Codon.CAG, Codon.GGA, Codon.CCC],  # QGP
            ),
            # Discontiguous CDS, minus strand, frame=2, with frameshift that leads to truncation
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_STRICT),
                    ),
                    [CDSFrame.TWO, CDSFrame.TWO],
                ),
                [Codon.CAG, Codon.GGA],  # QG
            ),
            # Discontiguous CDS, plus strand, with -1 bp programmed frameshift (overlapping interval)
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [9, 17],
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_STRICT),
                    ),
                    [CDSFrame.ZERO, CDSFrame.ONE],
                ),  # G gets repeated here
                [Codon.AGG, Codon.AAA, Codon.GGT, Codon.CCC, Codon.TGA],
            ),
        ],
    )
    def test_scan_codons(self, cds, expected):
        assert list(cds.scan_codons()) == expected

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # Contiguous CDS, plus strand, frame=0
            (
                CDSInterval(
                    SingleInterval(0, 9, Strand.PLUS, parent=Sequence("atacgatca", Alphabet.NT_STRICT)),
                    [CDSFrame.ZERO],
                ),
                [Codon.ATA, Codon.CGA, Codon.TCA],
            ),
            # Discontiguous CDS, plus strand, frame=1, codons don't reach end of CDS
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.PLUS,
                        parent=Sequence("aaacaaaagggacccaaaaaa", Alphabet.NT_STRICT),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                [Codon.CAG, Codon.GGA, Codon.CCC],  # QGP
            ),
        ],
    )
    def test_scan_codons_lower_case(self, cds, expected):
        assert list(cds.scan_codons()) == expected

    @pytest.mark.parametrize(
        "cds",
        [
            CDSInterval(
                SingleInterval(0, 9, Strand.PLUS, parent=Sequence("ANACGATCA", Alphabet.NT_EXTENDED_GAPPED)),
                [CDSFrame.ZERO],
            ),
            CDSInterval(
                CompoundInterval(
                    [2, 8],
                    [5, 17],
                    Strand.PLUS,
                    parent=Sequence("AANNNNAAGGGTACCCAAAAAA", Alphabet.NT_EXTENDED_GAPPED),
                ),
                [CDSFrame.ONE, CDSFrame.TWO],
            ),
        ],
    )
    def test_scan_codons_exception(self, cds):
        with pytest.raises(ValueError):
            _ = list(cds.scan_codons())

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # Contiguous CDS, plus strand, frame=0
            (
                CDSInterval(
                    SingleInterval(
                        0,
                        9,
                        Strand.PLUS,
                        parent=Sequence("ATACGATCA", Alphabet.NT_STRICT),
                    ),
                    [CDSFrame.ZERO],
                ),
                "IRS",
            ),
            # Discontiguous CDS, plus strand, frame=1, codons don't reach end of CDS
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", Alphabet.NT_STRICT),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                "QGP",
            ),
            # Discontiguous CDS, minus strand, frame=2, codons don't reach end of CDS
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_STRICT),
                    ),
                    [CDSFrame.TWO, CDSFrame.TWO],
                ),
                "QG",
            ),
        ],
    )
    def test_translate(self, cds, expected):
        assert str(cds.translate()) == expected

    def test_accessors(self):
        cds = CDSInterval(SingleInterval(0, 10, Strand.PLUS), [CDSFrame.ZERO])
        assert cds.start == cds.location.start
        assert cds.end == cds.location.end
        assert cds.strand == cds.location.strand

    @pytest.mark.parametrize(
        "cds,expected",
        [
            (
                CDSInterval(
                    SingleInterval(0, 9, Strand.PLUS, parent=Sequence("ATACGATGA", Alphabet.NT_STRICT)),
                    [CDSFrame.ZERO],
                ),
                True,
            ),
            (
                CDSInterval(
                    SingleInterval(0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", Alphabet.NT_STRICT)),
                    [CDSFrame.ZERO],
                ),
                False,
            ),
        ],
    )
    def test_has_valid_stop(self, cds, expected):
        assert cds.has_valid_stop == expected
