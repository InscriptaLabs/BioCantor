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
            (CDSFrame.NONE, 1, CDSFrame.NONE),
        ],
    )
    def test_shift(self, frame, shift, expected):
        assert frame.shift(shift) == expected


class TestCDS:

    alphabet = Alphabet.NT_STRICT

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
            # 2bp CDS
            (
                CDSInterval(
                    SingleInterval(0, 2, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet)),
                    [CDSFrame.ZERO],
                ),
                0,
            ),
            # Contiguous CDS, plus strand, frame=0
            (
                CDSInterval(
                    SingleInterval(0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet)),
                    [CDSFrame.ZERO],
                ),
                3,
            ),
            # Discontiguous CDS, plus strand, frame=1, codons don't reach end of CDS
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                3,
            ),
            # Discontiguous CDS, plus strand, with -1 bp programmed frameshift (overlapping interval)
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [9, 17],
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet),
                    ),
                    [CDSFrame.ZERO, CDSFrame.ONE],
                ),
                5,
            ),
        ],
    )
    def test_num_codons(self, cds, expected):
        assert cds.num_codons() == expected

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # Contiguous CDS, plus strand, frame=0
            (
                CDSInterval(
                    SingleInterval(0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet)),
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
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet),
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
                        parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet),
                    ),
                    [CDSFrame.ONE, CDSFrame.ZERO],
                ),
                [Codon.GGA, Codon.CCC],  # GP
            ),
            # Discontiguous CDS, plus strand, frame=1,
            # 1bp insertion inside exon 2 relative to some canonical genome and we want to maintain original frame
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8, 12],
                        [5, 11, 18],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet),
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
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet),
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
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet),
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
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet),
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
                    SingleInterval(0, 9, Strand.PLUS, parent=Sequence("atacgatca", alphabet)),
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
                        parent=Sequence("aaacaaaagggacccaaaaaa", alphabet),
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
            # 2bp CDS
            (
                CDSInterval(
                    SingleInterval(0, 2, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet)),
                    [CDSFrame.ZERO],
                ),
                [],
            ),
            # Contiguous CDS, plus strand, frame=0
            (
                CDSInterval(
                    SingleInterval(0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet)),
                    [CDSFrame.ZERO],
                ),
                [
                    SingleInterval(0, 3, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet)),  # ATA
                    SingleInterval(3, 6, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet)),  # CGA
                    SingleInterval(6, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet)),  # TCA
                ],
            ),
            # Discontiguous CDS, plus strand, frame=1, codons don't reach end of CDS
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                [
                    CompoundInterval(
                        [3, 8], [5, 9], Strand.PLUS, parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet)
                    ),  # CAG
                    SingleInterval(9, 12, Strand.PLUS, parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet)),  # GGA
                    SingleInterval(12, 15, Strand.PLUS, parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet)),  # CCC
                ],
            ),
            # Discontiguous CDS, plus strand, frame=1, 1bp deletion at start of exon 2
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [5, 16],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet),
                    ),
                    [CDSFrame.ONE, CDSFrame.ZERO],
                ),
                [
                    SingleInterval(8, 11, Strand.PLUS, parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet)),  # GGA
                    SingleInterval(11, 14, Strand.PLUS, parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet)),  # CCC
                ],
            ),
            # Discontiguous CDS, plus strand, frame=1,
            # 1bp insertion inside exon 2 relative to some canonical genome and we want to maintain original frame
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8, 12],
                        [5, 11, 18],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO, CDSFrame.TWO],
                ),
                [
                    CompoundInterval(
                        [3, 8], [5, 9], Strand.PLUS, parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet)
                    ),  # CAG
                    CompoundInterval(
                        [9, 12], [11, 13], Strand.PLUS, parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet)
                    ),  # GGA
                    SingleInterval(13, 16, Strand.PLUS, parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet)),  # CCC
                ],
            ),
            # Discontiguous CDS, minus strand, frame=2
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                [
                    SingleInterval(12, 15, Strand.MINUS, parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet)),  # CAG
                    SingleInterval(9, 12, Strand.MINUS, parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet)),  # GGA
                    CompoundInterval(
                        [3, 8], [5, 9], Strand.MINUS, parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet)
                    ),  # CCC
                ],
            ),
            # Discontiguous CDS, minus strand, frame=2, with frameshift that leads to truncation
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet),
                    ),
                    [CDSFrame.TWO, CDSFrame.TWO],
                ),
                [
                    SingleInterval(12, 15, Strand.MINUS, parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet)),  # CAG
                    SingleInterval(9, 12, Strand.MINUS, parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet)),  # GGA
                ],
            ),
            # Discontiguous CDS, plus strand, with -1 bp programmed frameshift (overlapping interval)
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 8],
                        [9, 17],
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet),
                    ),
                    [CDSFrame.ZERO, CDSFrame.ONE],
                ),
                [
                    SingleInterval(2, 5, Strand.PLUS, parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet)),  # AGG
                    SingleInterval(5, 8, Strand.PLUS, parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet)),  # AAA
                    CompoundInterval(
                        [8, 8], [9, 10], Strand.PLUS, parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet)
                    ),  # GGT, G gets repeated
                    SingleInterval(10, 13, Strand.PLUS, parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet)),  # CCC
                    SingleInterval(13, 16, Strand.PLUS, parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet)),  # TGA
                ],
            ),
            # Discontiguous CDS, plus strand, with +1 programmed frameshift that skips over a 1nt exon
            (
                CDSInterval(
                    CompoundInterval(
                        [2, 6, 8],
                        [5, 7, 16],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet),
                    ),
                    [CDSFrame.ZERO, CDSFrame.ZERO, CDSFrame.ZERO],
                ),
                [
                    SingleInterval(2, 5, Strand.PLUS, parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet)),  # GGA
                    SingleInterval(8, 11, Strand.PLUS, parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet)),  # GGA
                    SingleInterval(11, 14, Strand.PLUS, parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet)),  # CCC
                ],
            ),
        ],
    )
    def test_scan_codon_locations(self, cds, expected):
        assert list(cds.scan_codon_locations()) == expected

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
                        parent=Sequence("ATACGATCA", alphabet),
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
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet),
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
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet),
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
                    SingleInterval(0, 9, Strand.PLUS, parent=Sequence("ATACGATGA", alphabet)),
                    [CDSFrame.ZERO],
                ),
                True,
            ),
            (
                CDSInterval(
                    SingleInterval(0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet)),
                    [CDSFrame.ZERO],
                ),
                False,
            ),
        ],
    )
    def test_has_valid_stop(self, cds, expected):
        assert cds.has_valid_stop == expected
